"""
Iterative deformation engine: depth 0 (seeds) -> depth 1 -> ... -> depth 5.

Full pipeline at each depth:
1. Orbit decomposition (fields grouped by residual flavor symmetry)
2. Orbit-aware PE (Hilbert series with orbit fugacities)
3. Monomial expansion (PE degrees -> explicit monomials)
4. F-term reduction (GroebnerBasis via Mathematica)
5. Classify: relevant (R<2), super-relevant (R<4/3), unitarity-violating (R<=2/3)
6. Generate deformations (direct + flip) from each relevant/super-relevant operator
7. Deduplicate via graph isomorphism
8. Run a-maximization on candidates (parallel)
9. Repeat
"""

import json
import subprocess
import sys
import os
import time
import ast
import multiprocessing as mp
from fractions import Fraction

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from hilbert.hilbert import compute_hilbert_series
from hilbert.orbits import compute_orbits, orbits_to_pe_input
from hilbert.invariants import get_contraction_type
from hilbert.monomials import generate_all_operators
from deform.equivalence import deduplicate_theories

AMAX_PATH = "src/amax/FindCharges.wl"
REDUCE_PATH = "src/amax/ReduceOperators.wl"


def _amax_worker(args):
    """Worker function for parallel a-maximization."""
    dimG_, Tadj_, repinfo_, w_ = args
    return run_findcharges(dimG_, Tadj_, repinfo_, w_)


GROUP_DATA = {
    'A1': (3, 2), 'A2': (8, 3), 'A3': (15, 4), 'A4': (24, 5),
    'B2': (10, 3), 'B3': (21, 5),
    'C2': (10, 3), 'C3': (21, 4),
    'D4': (28, 6),
    'G2': (14, 4), 'F4': (52, 9),
}


def run_findcharges(dimG, Tadj, repinfo, w):
    """Call FindCharges via wolframscript. Returns parsed result dict."""
    ri_str = '{' + ', '.join(
        f'{{"{name}", {T_R}, {dim_r}, {n}}}'
        for name, T_R, dim_r, n in repinfo
    ) + '}'
    w_str = '{' + ', '.join(f'"{t}"' for t in w) + '}'

    mcode = f"""
Get["{AMAX_PATH}"];
result = FindCharges[{dimG}, {Tadj}, {ri_str}, {w_str}];
Print[result["consistency"]];
Print[N[result["a"], 15]];
Print[N[result["c"], 15]];
Do[
  Print[k, " -> ", N[result[k], 15]],
  {{k, Select[Keys[result],
    StringQ[#] && !MemberQ[{{"n","w","consistency","method","a","c",
      "rational","computation_time","global"}}, #] &]}}
];
"""
    try:
        proc = subprocess.run(
            ['wolframscript', '-code', mcode],
            capture_output=True, encoding='UTF-8', timeout=300
        )
    except subprocess.TimeoutExpired:
        return {'consistency': 'timeout'}

    lines = [l.strip() for l in proc.stdout.split('\n') if l.strip()]
    if len(lines) < 3:
        return {'consistency': 'error'}

    result = {'consistency': lines[0], 'a': lines[1], 'c': lines[2]}
    for line in lines[3:]:
        if ' -> ' in line:
            k, v = line.split(' -> ', 1)
            result[k.strip()] = v.strip()
    return result


def run_reduce_operators(monomial_strs, field_names, w_terms, r_charges):
    """Call ReduceAndClassify via wolframscript. Returns classified operators."""
    mono_mma = '{' + ','.join(f'"{m}"' for m in monomial_strs) + '}'
    fields_mma = '{' + ','.join(f'"{f}"' for f in field_names) + '}'
    w_mma = '{' + ','.join(f'"{t}"' for t in w_terms) + '}'
    rc_mma = '<|' + ','.join(f'"{k}"->{v}' for k, v in r_charges.items()) + '|>'

    mcode = f"""
Get["{REDUCE_PATH}"];
result = ReduceAndClassify[{mono_mma}, {fields_mma}, {w_mma}, {rc_mma}];
Print[result["n_relevant"]];
Print[result["n_flippable"]];
Print[result["n_unitarity_violating"]];
Do[Print[r[[1]], " ||| ", r[[2]]], {{r, result["relevant_with_R"]}}];
"""
    try:
        proc = subprocess.run(
            ['wolframscript', '-code', mcode],
            capture_output=True, encoding='UTF-8', timeout=120
        )
    except subprocess.TimeoutExpired:
        return None

    lines = [l.strip() for l in proc.stdout.split('\n')
             if l.strip() and not l.startswith('Failed')]
    if len(lines) < 3:
        return None

    n_relevant = int(lines[0])
    n_flippable = int(lines[1])
    n_unitarity = int(lines[2])

    relevant_ops = []
    flippable_ops = []
    for line in lines[3:]:
        if ' ||| ' in line:
            mono, r_val = line.split(' ||| ')
            R = float(r_val.split('`')[0]) if '`' in r_val else float(r_val)
            relevant_ops.append({'monomial': mono, 'R': R})
            if 2/3 + 1e-10 < R < 4/3:
                flippable_ops.append({'monomial': mono, 'R': R})

    return {
        'n_relevant': n_relevant,
        'n_flippable': n_flippable,
        'n_unitarity_violating': n_unitarity,
        'relevant': relevant_ops,
        'flippable': flippable_ops,
    }


def parse_float(s):
    """Parse Mathematica-style float (with backtick precision)."""
    return float(s.split('`')[0])


def get_theory_operators(group, theory, all_reps, hs_order=6):
    """Full operator enumeration for a theory using orbit-aware PE.

    Returns:
        dict with 'relevant', 'flippable', 'has_unitarity_violation', or None on failure
    """
    matter = theory['matter']
    w = theory.get('w', [])
    rcharges_raw = theory.get('rcharges', {})
    stored_repinfo = theory.get('repinfo', [])

    # Parse R-charges and build field lists using repinfo names
    r_charges = {}
    all_fields = []
    rep_type_to_dynkin = {}
    singlet_prefixes = set()

    # Build from stored repinfo if available (preserves naming like M0, r0, etc.)
    if stored_repinfo:
        for ri_idx, (name, T_R, dim_r, n) in enumerate(stored_repinfo):
            is_singlet = (str(T_R) == '0' and dim_r == 1)
            if is_singlet:
                singlet_prefixes.add(name)
            for copy in range(1, n + 1):
                fname = f"{name}{copy}"
                all_fields.append(fname)
                rk = rcharges_raw.get(fname)
                r_charges[fname] = parse_float(rk) if rk else 0.5
            # Find Dynkin label from matter
            # Match repinfo to matter: repinfo may have more entries (R + Rbar for complex)
            dl = None
            for midx, (label, mn, reality, md) in enumerate(matter):
                if f"r{midx}" == name or name == f"r{midx}b":
                    dl = label
                    break
                # Also check singlet names like M0, M1, etc.
                if is_singlet and all(x == 0 for x in label):
                    dl = label
                    break
            if dl is None:
                dl = [0]  # fallback: singlet
            rep_type_to_dynkin[name] = dl
    else:
        # Build from matter (depth 0 seeds)
        for midx, (label, n, reality, dim_r) in enumerate(matter):
            rep_type = f"r{midx}"
            T_R_data = None
            for r in all_reps[group]['reps']:
                if tuple(r['dynkin_label']) == tuple(label):
                    T_R_data = r['T_R']
                    break
            if T_R_data is None and all(x == 0 for x in label):
                T_R_data = '0'
                singlet_prefixes.add(rep_type)

            if reality == 'complex':
                for copy in range(1, n + 1):
                    fname = f"{rep_type}{copy}"
                    all_fields.append(fname)
                    rk = rcharges_raw.get(fname)
                    r_charges[fname] = parse_float(rk) if rk else 0.5
                rep_type_to_dynkin[rep_type] = label
                rep_type_b = f"r{midx}b"
                conj_label = label
                for r in all_reps[group]['reps']:
                    if tuple(r['dynkin_label']) == tuple(label) and r.get('contragr'):
                        conj_label = r['contragr']
                        break
                for copy in range(1, n + 1):
                    fname = f"{rep_type_b}{copy}"
                    all_fields.append(fname)
                    rk = rcharges_raw.get(fname)
                    r_charges[fname] = parse_float(rk) if rk else 0.5
                rep_type_to_dynkin[rep_type_b] = conj_label
            else:
                is_singlet = all(x == 0 for x in label)
                if is_singlet:
                    singlet_prefixes.add(rep_type)
                for copy in range(1, n + 1):
                    fname = f"{rep_type}{copy}"
                    all_fields.append(fname)
                    rk = rcharges_raw.get(fname)
                    r_charges[fname] = parse_float(rk) if rk else 0.5
                rep_type_to_dynkin[rep_type] = label

    # 1. Orbit decomposition
    orbits = compute_orbits(all_fields, w, singlet_prefixes=singlet_prefixes)

    # 2. Orbit-aware PE
    pe_labels, pe_mults, pe_oidx = orbits_to_pe_input(orbits, rep_type_to_dynkin)

    if not pe_labels:
        # All singlets — no gauge-charged fields
        pe_degrees = []
    else:
        hs = compute_hilbert_series(group, pe_labels, pe_mults, order=hs_order, compute_pl=False)
        pe_degrees = []
        for item_str in hs['pe_raw']:
            matrix = ast.literal_eval(item_str)
            rows = matrix if isinstance(matrix[0], list) else [matrix]
            for row in rows:
                pe_degrees.append((tuple(row[:-1]), row[-1]))

    # 3. Build orbit info for monomial expansion
    gauge_orbits = []
    for pi, oi in enumerate(pe_oidx):
        o = orbits[oi]
        dl = rep_type_to_dynkin.get(o['rep_type'], [0])
        ct = get_contraction_type(group, dl) if not o['is_singlet'] else 'any'
        gauge_orbits.append({
            'fields': o['fields'],
            'rep_idx': pi,
            'contraction': ct,
            'is_singlet': False,
        })

    singlet_fields = [f for o in orbits if o['is_singlet'] for f in o['fields']]

    # 4. Generate monomials (before F-term)
    ops = generate_all_operators(gauge_orbits, pe_degrees, {}, r_charges,
                                 singlet_fields=singlet_fields, max_R=2.0)

    if not ops and not singlet_fields:
        return {'relevant': [], 'flippable': [], 'has_unitarity_violation': False}

    # 5. F-term reduction via Mathematica
    mono_strs = [op['monomial'] for op in ops]
    if not mono_strs:
        return {'relevant': [], 'flippable': [], 'has_unitarity_violation': False}

    reduced = run_reduce_operators(mono_strs, all_fields, w, r_charges)
    if reduced is None:
        return None

    has_viol = reduced['n_unitarity_violating'] > 0

    return {
        'relevant': reduced['relevant'],
        'flippable': reduced['flippable'],
        'has_unitarity_violation': has_viol,
        'n_relevant': reduced['n_relevant'],
        'n_flippable': reduced['n_flippable'],
    }


def iterate_depth(group, seeds_at_depth, all_reps, depth, max_depth=5, hs_order=6):
    """Process one depth level: generate candidates at depth+1."""
    if depth >= max_depth:
        return []

    dimG, Tadj = GROUP_DATA[group]
    candidates = []

    for theory in seeds_at_depth:
        if theory.get('consistency', 'consistent') != 'consistent':
            continue

        # Get operators via orbit-aware PE + F-term reduction
        op_result = get_theory_operators(group, theory, all_reps, hs_order)
        if op_result is None:
            continue

        if op_result['has_unitarity_violation']:
            theory['consistency'] = 'operator decoupled'
            continue

        matter = theory['matter']
        w = theory.get('w', [])
        repinfo = theory.get('repinfo', [])

        # If repinfo not stored, rebuild it
        if not repinfo:
            for midx, (label, n, reality, dim_r) in enumerate(matter):
                T_R_data = '0'
                for r in all_reps[group]['reps']:
                    if tuple(r['dynkin_label']) == tuple(label):
                        T_R_data = r['T_R']
                        break
                rep_type = f"r{midx}"
                if reality == 'complex':
                    conj_label = label
                    for r in all_reps[group]['reps']:
                        if tuple(r['dynkin_label']) == tuple(label) and r.get('contragr'):
                            conj_label = r['contragr']
                            break
                    repinfo.append((rep_type, str(Fraction(T_R_data)), dim_r, n))
                    repinfo.append((f"r{midx}b", str(Fraction(T_R_data)), dim_r, n))
                else:
                    repinfo.append((rep_type, str(Fraction(T_R_data)), dim_r, n))

        # Generate direct deformations
        for op in op_result['relevant']:
            new_w = w + [op['monomial']]
            candidates.append({
                'repinfo': list(repinfo),
                'matter': list(matter),
                'w': new_w,
                'parent_desc': theory.get('description', ''),
                'deform_type': 'direct',
                'deform_op': op['monomial'],
                'deform_R': op['R'],
            })

        # Generate flip deformations
        for op in op_result['flippable']:
            flip_idx = len(matter)
            rank = len(matter[0][0]) if matter else 1
            singlet_label = [0] * rank
            singlet_name = f"M{flip_idx}"
            new_repinfo = list(repinfo) + [(singlet_name, '0', 1, 1)]
            new_matter = list(matter) + [(singlet_label, 1, 'real', 1)]
            new_w = w + [f"{singlet_name}1*{op['monomial']}"]
            candidates.append({
                'repinfo': new_repinfo,
                'matter': new_matter,
                'w': new_w,
                'parent_desc': theory.get('description', ''),
                'deform_type': 'flip',
                'deform_op': op['monomial'],
                'deform_R': op['R'],
            })

    # Deduplicate
    n_before = len(candidates)
    candidates = deduplicate_theories(candidates)
    n_after = len(candidates)

    print(f"  Depth {depth} -> {depth+1}: {n_before} candidates, "
          f"{n_after} after dedup, from {len(seeds_at_depth)} theories")

    # Run a-maximization in parallel
    n_workers = min(mp.cpu_count(), len(candidates)) if candidates else 1

    work_items = [
        (dimG, Tadj, cand['repinfo'], cand['w'])
        for cand in candidates
    ]

    if work_items:
        print(f"  Running a-maximization on {len(candidates)} candidates "
              f"with {n_workers} workers...")
        with mp.Pool(n_workers) as pool:
            amax_results = pool.map(_amax_worker, work_items)
    else:
        amax_results = []

    # Collect consistent results
    results = []
    for cand, result in zip(candidates, amax_results):
        if result['consistency'] == 'consistent':
            result['matter'] = cand['matter']
            result['w'] = cand['w']
            result['repinfo'] = cand['repinfo']
            result['depth'] = depth + 1
            result['parent'] = cand['parent_desc']
            result['deform_type'] = cand['deform_type']
            result['deform_op'] = cand['deform_op']
            result['description'] = f"{cand['parent_desc']} + {cand['deform_type']}({cand['deform_op']})"
            result['rcharges'] = {k: v for k, v in result.items()
                                  if k not in ('consistency', 'a', 'c', 'matter', 'w',
                                              'repinfo', 'depth', 'parent', 'deform_type',
                                              'deform_op', 'description')}
            results.append(result)

    print(f"  {len(results)}/{len(candidates)} consistent")
    return results


if __name__ == '__main__':
    with open('data/all_reps.json') as f:
        all_reps = json.load(f)
    with open('data/consistent_seeds.json') as f:
        seeds_data = json.load(f)

    group = sys.argv[1] if len(sys.argv) > 1 else 'G2'
    max_depth = int(sys.argv[2]) if len(sys.argv) > 2 else 1
    seed_idx = int(sys.argv[3]) if len(sys.argv) > 3 else None

    if group not in seeds_data:
        print(f"No seeds for {group}")
        sys.exit(1)

    print(f"Iterating {group} to depth {max_depth}")

    all_seeds = seeds_data[group]
    if seed_idx is not None:
        all_seeds = [all_seeds[seed_idx]]
        print(f"  Single seed #{seed_idx}: {all_seeds[0]['description']}")

    current_depth_theories = []
    for seed in all_seeds:
        theory = dict(seed)
        theory['w'] = []
        theory['depth'] = 0
        theory['consistency'] = 'consistent'
        current_depth_theories.append(theory)

    all_theories = list(current_depth_theories)

    for d in range(max_depth):
        print(f"\n{'='*60}")
        print(f"  Processing depth {d} -> {d+1}")
        print(f"{'='*60}")

        next_theories = iterate_depth(
            group, current_depth_theories, all_reps, d,
            max_depth=max_depth, hs_order=6
        )

        print(f"  Depth {d+1}: {len(next_theories)} consistent theories")
        all_theories.extend(next_theories)
        current_depth_theories = next_theories

    print(f"\nTotal: {len(all_theories)} theories (depth 0 to {max_depth})")
    for d in range(max_depth + 1):
        count = sum(1 for t in all_theories if t.get('depth', 0) == d)
        print(f"  Depth {d}: {count}")
