"""
Iterative deformation engine: depth 0 (seeds) -> depth 1 -> ... -> depth 5.

For each consistent theory at depth n:
1. Enumerate relevant (R<2) and super-relevant (R<4/3) GIO degrees from PE
2. For each relevant degree: direct deformation dW = O (one per inequivalent degree)
3. For each super-relevant degree: flip dW = M*O (adds singlet M)
4. Run a-maximization on each candidate at depth n+1
5. Compute PE + classify operators for consistent results
6. Repeat until depth 5

For W=0 seeds, operators at the same degree are related by flavor symmetry,
so there is one inequivalent deformation per degree. For W!=0, deduplication
via graph isomorphism (toGlobalGraph) would be needed — simplified here.
"""

import json
import subprocess
import sys
import os
import time
from fractions import Fraction

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from hilbert.hilbert import compute_hilbert_series
from deform.equivalence import deduplicate_theories

AMAX_PATH = "src/amax/FindCharges.wl"

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


def parse_float(s):
    """Parse Mathematica-style float (with backtick precision)."""
    return float(s.split('`')[0])


def get_pe_operators(lie_group, dynkin_labels, multiplicities, order=6):
    """Compute PE and return parsed operator list."""
    import ast
    hs = compute_hilbert_series(lie_group, dynkin_labels, multiplicities, order, compute_pl=False)
    operators = []
    for item_str in hs['pe_raw']:
        matrix = ast.literal_eval(item_str)
        if isinstance(matrix[0], list):
            for row in matrix:
                operators.append((tuple(row[:-1]), row[-1]))
        else:
            operators.append((tuple(matrix[:-1]), matrix[-1]))
    return operators


def classify_pe_ops(operators, rcharges_per_type):
    """Classify PE operators by R-charge. Returns relevant and super-relevant degrees."""
    relevant = []       # R < 2
    super_relevant = []  # R < 4/3

    for degree, mult in operators:
        if all(d == 0 for d in degree):
            continue
        R = sum(d * r for d, r in zip(degree, rcharges_per_type))
        if R < 2:
            relevant.append({'degree': degree, 'mult': mult, 'R': R})
        if R < 4/3:
            super_relevant.append({'degree': degree, 'mult': mult, 'R': R})

    return relevant, super_relevant


def build_monomial_str(degree, rep_names, copy_start=1):
    """Build a representative monomial string from a degree vector.

    E.g. degree=(2,1) with names=['q','phi'] -> 'q1*q2*phi1'
    """
    parts = []
    for i, d in enumerate(degree):
        for j in range(d):
            parts.append(f"{rep_names[i]}{copy_start + j}")
    return '*'.join(parts)


def iterate_depth(group, seeds_at_depth, all_reps, depth, max_depth=5, hs_order=6):
    """Process one depth level: generate candidates at depth+1.

    Args:
        group: LiE group name
        seeds_at_depth: list of theory dicts at current depth
        all_reps: rep data from all_reps.json
        depth: current depth
        max_depth: stop at this depth
        hs_order: Hilbert series truncation order

    Returns:
        list of consistent theory dicts at depth+1
    """
    if depth >= max_depth:
        return []

    dimG, Tadj = GROUP_DATA[group]
    rep_data = {tuple(r['dynkin_label']): r for r in all_reps[group]['reps']}

    candidates = []  # (repinfo, w_list, description)

    for theory in seeds_at_depth:
        if theory.get('consistency', 'consistent') != 'consistent':
            continue

        matter = theory['matter']
        w = theory.get('w', [])
        rcharges = theory.get('rcharges', {})

        # Build rep info for Hilbert series (gauge-charged only)
        hs_labels = []
        hs_mults = []
        hs_names = []
        hs_rcharges = []

        repinfo = []  # for FindCharges
        rep_names_for_monomial = []

        for midx, (label, n, reality, dim_r) in enumerate(matter):
            r_key = f"r{midx}1"
            r_charge = 0.5
            for k, v in rcharges.items():
                if k == r_key or (k.startswith(f"r{midx}") and 'b' not in k):
                    try:
                        r_charge = parse_float(v)
                    except:
                        pass
                    break

            T_R_str = rep_data.get(tuple(label), {}).get('T_R', '1/2')
            T_R = str(Fraction(T_R_str))

            if reality == 'complex':
                name_r = f"r{midx}"
                name_rb = f"r{midx}b"
                repinfo.append((name_r, T_R, dim_r, n))
                repinfo.append((name_rb, T_R, dim_r, n))

                conj = rep_data.get(tuple(label), {}).get('contragr', label)
                hs_labels.append(label)
                hs_mults.append(n)
                hs_names.append(name_r)
                hs_rcharges.append(r_charge)
                hs_labels.append(conj if conj else label)
                hs_mults.append(n)
                hs_names.append(name_rb)
                hs_rcharges.append(r_charge)
                rep_names_for_monomial.extend([name_r, name_rb])
            else:
                name_r = f"r{midx}"
                repinfo.append((name_r, T_R, dim_r, n))
                hs_labels.append(label)
                hs_mults.append(n)
                hs_names.append(name_r)
                hs_rcharges.append(r_charge)
                rep_names_for_monomial.append(name_r)

        # Compute PE and classify operators
        operators = get_pe_operators(group, hs_labels, hs_mults, order=hs_order)
        relevant, super_relevant = classify_pe_ops(operators, hs_rcharges)

        # Check for unitarity violations
        has_unitarity_violation = any(
            sum(d * r for d, r in zip(op['degree'], hs_rcharges)) < 2/3
            for op in relevant
        )
        if has_unitarity_violation:
            theory['consistency'] = 'operator decoupled'
            continue

        # Generate direct deformations: one per relevant degree
        for op in relevant:
            monomial = build_monomial_str(op['degree'], hs_names)
            new_w = w + [monomial]
            candidates.append({
                'repinfo': list(repinfo),
                'matter': list(matter),
                'w': new_w,
                'parent_desc': theory.get('description', ''),
                'deform_type': 'direct',
                'deform_op': monomial,
                'deform_R': op['R'],
            })

        # Generate flip deformations: one per super-relevant degree
        n_singlets = sum(1 for _, _, r, _ in matter if r == 'real' and _ == 1)  # rough count
        flip_idx = len(matter)  # index for new singlet

        for op in super_relevant:
            monomial = build_monomial_str(op['degree'], hs_names)
            singlet_name = f"M{flip_idx}"
            flip_term = f"{singlet_name}1*{monomial}"

            # New matter: add one singlet
            new_matter = list(matter) + [([0] * len(matter[0][0]) if matter else [0], 1, 'real', 1)]
            new_repinfo = list(repinfo) + [(singlet_name, 0, 1, 1)]
            new_w = w + [flip_term]

            candidates.append({
                'repinfo': new_repinfo,
                'matter': new_matter,
                'w': new_w,
                'parent_desc': theory.get('description', ''),
                'deform_type': 'flip',
                'deform_op': monomial,
                'deform_R': op['R'],
            })

    # Deduplicate candidates before running a-maximization
    n_before = len(candidates)
    candidates = deduplicate_theories(candidates)
    n_after = len(candidates)

    print(f"  Depth {depth} -> {depth+1}: {n_before} candidates, "
          f"{n_after} after dedup, from {len(seeds_at_depth)} theories")

    # Run a-maximization on each candidate
    results = []
    for i, cand in enumerate(candidates):
        result = run_findcharges(dimG, Tadj, cand['repinfo'], cand['w'])

        if result['consistency'] == 'consistent':
            result['matter'] = cand['matter']
            result['w'] = cand['w']
            result['depth'] = depth + 1
            result['parent'] = cand['parent_desc']
            result['deform_type'] = cand['deform_type']
            result['deform_op'] = cand['deform_op']
            result['description'] = f"{cand['parent_desc']} + {cand['deform_type']}({cand['deform_op']})"
            result['rcharges'] = {k: v for k, v in result.items()
                                  if k not in ('consistency', 'a', 'c', 'matter', 'w',
                                              'depth', 'parent', 'deform_type', 'deform_op',
                                              'description')}
            results.append(result)

        if (i + 1) % 10 == 0 or i == len(candidates) - 1:
            print(f"    [{i+1}/{len(candidates)}] {len(results)} consistent so far")

    return results


if __name__ == '__main__':
    with open('data/all_reps.json') as f:
        all_reps = json.load(f)
    with open('data/consistent_seeds.json') as f:
        seeds_data = json.load(f)

    group = sys.argv[1] if len(sys.argv) > 1 else 'G2'
    max_depth = int(sys.argv[2]) if len(sys.argv) > 2 else 1

    if group not in seeds_data:
        print(f"No seeds for {group}")
        sys.exit(1)

    print(f"Iterating {group} to depth {max_depth}")
    dimG, Tadj = GROUP_DATA[group]

    # Convert seeds to the format expected by iterate_depth
    current_depth_theories = []
    for seed in seeds_data[group]:
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
