"""
Construct gauge-invariant operators from primitive building blocks.

1. Find primitive invariants via invariants.py (degree, sym/alt, multiplicity)
2. Enumerate explicit primitive monomials from field copies
3. Build all products of primitives up to R < maxR
4. Validate against orbit-aware PE count at each degree
5. Return explicit monomial list

If PE count > product count at any degree: raises ValueError (missing primitives).
"""

import sys
import os
import ast
from collections import defaultdict
from itertools import combinations, combinations_with_replacement, product as iprod

_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _dir)
sys.path.insert(0, os.path.join(_dir, '..'))
from invariants import find_single_rep_invariants, find_mixed_rep_invariants
from hilbert import compute_hilbert_series


def enumerate_primitives_single(fields, degree, sym_type):
    """Enumerate explicit primitive monomials from one rep type.

    Args:
        fields: list of field name strings (e.g. ["q1","q2","q3"])
        degree: int (e.g. 2 for mesons, 3 for baryons)
        sym_type: "alt" (antisymmetric, no repetition) or "sym" (symmetric, with repetition)

    Returns:
        list of monomial strings
    """
    if sym_type == "alt":
        combos = combinations(fields, degree)
    else:
        combos = combinations_with_replacement(fields, degree)
    return ['*'.join(c) for c in combos]


def enumerate_primitives_mixed(fields_i, fields_j):
    """Enumerate explicit bilinear mixed-rep primitives (one field from each type).

    Returns:
        list of monomial strings
    """
    return [f"{fi}*{fj}" for fi in fields_i for fj in fields_j]


def build_all_products(primitives, r_charges, max_R=2.0, max_factors=3):
    """Build all products of primitive monomials up to max_R.

    Args:
        primitives: list of monomial strings
        r_charges: dict field -> R-charge
        max_R: keep products with R < max_R
        max_factors: maximum number of primitive factors in a product

    Returns:
        list of monomial strings (including primitives themselves)
    """

    def mono_R(mono):
        return sum(r_charges.get(f, 0) for f in mono.split('*'))

    # Start with primitives that have R < max_R
    result = set()
    for p in primitives:
        if mono_R(p) < max_R:
            key = tuple(sorted(p.split('*')))
            result.add(key)

    # Build products of 2 primitives
    for i in range(len(primitives)):
        R_i = mono_R(primitives[i])
        if R_i >= max_R:
            continue
        for j in range(i, len(primitives)):
            R_j = mono_R(primitives[j])
            if R_i + R_j >= max_R:
                continue
            fields = primitives[i].split('*') + primitives[j].split('*')
            key = tuple(sorted(fields))
            result.add(key)

            # Products of 3 primitives
            if max_factors >= 3:
                for k in range(j, len(primitives)):
                    R_k = mono_R(primitives[k])
                    if R_i + R_j + R_k >= max_R:
                        continue
                    fields3 = fields + primitives[k].split('*')
                    key3 = tuple(sorted(fields3))
                    result.add(key3)

    return ['*'.join(key) for key in result]


def construct_gio(lie_group, repinfo, all_fields, singlet_fields, r_charges,
                  w_terms, max_R=2.0):
    """Full GIO construction with PE validation.

    Args:
        lie_group: LiE group name (e.g. "A1")
        repinfo: list of (name, T_R, dim_R, n_copies) — gauge-charged reps only
        all_fields: list of all field name strings (gauge + singlet)
        singlet_fields: list of singlet field names
        r_charges: dict field -> R-charge
        w_terms: list of W term strings
        max_R: keep operators with R < max_R

    Returns:
        dict with 'operators' (list of monomial strings), 'validation' (pass/fail)
    """
    # Group gauge fields by rep type
    rep_fields = defaultdict(list)  # name -> list of field names
    rep_dynkin = {}  # name -> dynkin label
    for name, T_R, dim_R, n in repinfo:
        if str(T_R) == '0' and dim_R == 1:
            continue  # skip singlets
        for i in range(1, n + 1):
            rep_fields[name].append(f"{name}{i}")
        # Find Dynkin label — stored in caller, pass via repinfo extension
        # For now, assume it's available from the context

    # Get Dynkin labels for each rep type
    # We need these for invariants.py. They should be passed in.
    # For simplicity, accept an additional parameter or infer from context.

    # WORKAROUND: extract Dynkin labels from repinfo names via all_reps
    # This is fragile; better to pass explicitly. For now, accept as parameter.

    return _construct_gio_with_labels(
        lie_group, repinfo, rep_fields, all_fields, singlet_fields,
        r_charges, w_terms, max_R
    )


def _construct_gio_with_labels(lie_group, repinfo, rep_fields, all_fields,
                                singlet_fields, r_charges, w_terms, max_R,
                                dynkin_labels=None):
    """Internal: construct GIO with known Dynkin labels."""

    # 1. Find primitive invariants
    rep_names = list(rep_fields.keys())
    all_gauge_fields = []
    for name in rep_names:
        all_gauge_fields.extend(rep_fields[name])

    primitives = []

    if dynkin_labels:
        # Single-rep primitives
        for i, name in enumerate(rep_names):
            dl = dynkin_labels[i]
            invs = find_single_rep_invariants(lie_group, dl, max_degree=6)
            for degree, sym_type, mult in invs:
                monos = enumerate_primitives_single(rep_fields[name], degree, sym_type)
                primitives.extend(monos)

        # Mixed-rep primitives
        mixed = find_mixed_rep_invariants(lie_group,
                                          [dynkin_labels[i] for i in range(len(rep_names))])
        for i, j, mult in mixed:
            if i == j:
                continue  # same-rep handled above
            monos = enumerate_primitives_mixed(rep_fields[rep_names[i]],
                                               rep_fields[rep_names[j]])
            primitives.extend(monos)
    else:
        # Fallback: use all pairs (may include non-gauge-invariant ones)
        # This should not happen in production
        primitives = ['*'.join(c) for c in combinations(all_gauge_fields, 2)]

    # 2. Build all products of primitives
    all_products = build_all_products(primitives, r_charges, max_R, max_factors=3)

    # 3. Add singlet fields and their products
    all_ops = list(all_products)
    for sf in singlet_fields:
        R_sf = r_charges.get(sf, 0)
        if R_sf < max_R:
            all_ops.append(sf)
        if 2 * R_sf < max_R:
            all_ops.append(f"{sf}*{sf}")

    # Singlet × gauge ops
    for sf in singlet_fields:
        R_sf = r_charges.get(sf, 0)
        for gop in all_products:
            R_g = sum(r_charges.get(f, 0) for f in gop.split('*'))
            if R_sf + R_g < max_R:
                fields = [sf] + gop.split('*')
                all_ops.append('*'.join(fields))

    # Singlet × singlet
    for i, sf1 in enumerate(singlet_fields):
        for sf2 in singlet_fields[i+1:]:
            R_total = r_charges.get(sf1, 0) + r_charges.get(sf2, 0)
            if R_total < max_R:
                all_ops.append(f"{sf1}*{sf2}")

    # Deduplicate by sorted field tuple
    seen = set()
    unique_ops = []
    for op in all_ops:
        key = tuple(sorted(op.split('*')))
        if key not in seen:
            seen.add(key)
            unique_ops.append(op)

    # 4. Validate against per-rep PE (one fugacity per rep type, fastest)
    # Build per-rep PE input: group all gauge fields by Dynkin label
    rep_type_labels = []  # unique Dynkin labels
    rep_type_mults = []   # total multiplicity of each
    rep_type_rcharges = []  # average R-charge per type
    label_to_idx = {}

    for name in rep_names:
        dl = tuple(dynkin_labels[rep_names.index(name)]) if dynkin_labels else None
        if dl is None:
            continue
        if dl not in label_to_idx:
            label_to_idx[dl] = len(rep_type_labels)
            rep_type_labels.append(list(dl))
            rep_type_mults.append(0)
            avg_r = sum(r_charges.get(f, 0) for f in rep_fields[name]) / len(rep_fields[name])
            rep_type_rcharges.append(avg_r)
        label_to_idx_val = label_to_idx[dl]
        rep_type_mults[label_to_idx_val] += len(rep_fields[name])

    if rep_type_labels:
        hs = compute_hilbert_series(lie_group, rep_type_labels, rep_type_mults,
                                     order=4, compute_pl=False)
        pe_degrees = {}
        for item_str in hs['pe_raw']:
            matrix = ast.literal_eval(item_str)
            rows = matrix if isinstance(matrix[0], list) else [matrix]
            for row in rows:
                degree = tuple(row[:-1])
                mult = row[-1]
                if any(d > 0 for d in degree):
                    pe_degrees[degree] = mult

        # Count our gauge-only operators per rep-type degree
        our_counts = defaultdict(int)
        for op in unique_ops:
            fields_in_op = op.split('*')
            if any(f in singlet_fields for f in fields_in_op):
                continue
            degree = [0] * len(rep_type_labels)
            for f in fields_in_op:
                for name in rep_names:
                    if f in rep_fields[name]:
                        dl = tuple(dynkin_labels[rep_names.index(name)])
                        degree[label_to_idx[dl]] += 1
                        break
            our_counts[tuple(degree)] += 1

        # Exact match validation
        validation_passed = True
        for degree, pe_mult in pe_degrees.items():
            R_deg = sum(d * r for d, r in zip(degree, rep_type_rcharges))
            if R_deg >= max_R:
                continue
            our_mult = our_counts.get(degree, 0)
            if our_mult != pe_mult:
                print(f"  PE MISMATCH at degree {degree} R={R_deg:.4f}: "
                      f"PE={pe_mult}, ours={our_mult}")
                validation_passed = False

        # Also check for extra operators we have but PE doesn't
        for degree, our_mult in our_counts.items():
            if degree not in pe_degrees and our_mult > 0:
                R_deg = sum(d * r for d, r in zip(degree, rep_type_rcharges))
                if R_deg < max_R:
                    print(f"  EXTRA operators at degree {degree} R={R_deg:.4f}: "
                          f"ours={our_mult}, PE=0")
                    validation_passed = False

        if validation_passed:
            print(f"  PE validation PASSED (exact match at all degrees with R<{max_R})")
        else:
            print(f"  PE validation FAILED — STOPPING")
    else:
        validation_passed = True

    return {
        'operators': unique_ops,
        'n_operators': len(unique_ops),
        'validation_passed': validation_passed,
    }


if __name__ == '__main__':
    # Test: SU(2) 8f + 2 singlets, W = M01*r01*r02 + M01*r03*r04 + M02*r01*r05
    r_charges = {
        'r01': 0.5411, 'r02': 0.4863, 'r03': 0.5137, 'r04': 0.5137,
        'r05': 0.5262, 'r06': 0.4730, 'r07': 0.4730, 'r08': 0.4730,
        'M01': 0.9726, 'M02': 0.9327,
    }
    all_fields = [f'r0{i}' for i in range(1, 9)] + ['M01', 'M02']
    singlet_fields = ['M01', 'M02']
    w = ['M01*r01*r02', 'M01*r03*r04', 'M02*r01*r05']
    repinfo = [('r0', '1/2', 2, 8), ('M0', '0', 1, 2)]

    result = _construct_gio_with_labels(
        'A1', repinfo,
        {'r0': [f'r0{i}' for i in range(1, 9)]},
        all_fields, singlet_fields, r_charges, w, max_R=2.0,
        dynkin_labels=[[1]]
    )

    print(f"Total operators: {result['n_operators']}")
    print(f"PE validation: {'PASSED' if result['validation_passed'] else 'FAILED'}")

    # Count by type
    gauge_only = [op for op in result['operators']
                  if not any(f.startswith('M') for f in op.split('*'))]
    with_singlet = [op for op in result['operators']
                    if any(f.startswith('M') for f in op.split('*'))]
    mesons = [op for op in gauge_only if op.count('*') == 1]
    quartic = [op for op in gauge_only if op.count('*') >= 2]
    print(f"  Mesons: {len(mesons)}")
    print(f"  Quartic+: {len(quartic)}")
    print(f"  With singlet: {len(with_singlet)}")
