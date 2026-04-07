"""
Analyze a seed theory: a-maximization + Hilbert series + operator classification.

For each consistent seed from consistent_seeds.json:
1. Extract R-charges from a-maximization (already computed)
2. Compute Hilbert series via PE to enumerate GIOs
3. Assign R-charges to each GIO degree
4. Classify: relevant (R<2), super-relevant (R<4/3), unitarity-violating (R<2/3)
5. Output the full seed catalog with operator counts
"""

import json
import sys
import os
from fractions import Fraction

# Add parent to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from hilbert.hilbert import compute_hilbert_series


def parse_pe_raw(pe_raw_list):
    """Parse PE raw output from LiE into list of (degree_tuple, multiplicity).

    Each entry in pe_raw_list is a string like '[[0,2,1]]' or '[[1,3,1],[1,5,1]]'
    representing a LiE matrix. Each row is [n1, n2, ..., mult].
    """
    import ast
    if not pe_raw_list:
        return []
    results = []
    for item_str in pe_raw_list:
        matrix = ast.literal_eval(item_str) if isinstance(item_str, str) else item_str
        if isinstance(matrix, list):
            # Could be a single row [[a,b,c]] or multi-row [[a,b,c],[d,e,f]]
            if matrix and isinstance(matrix[0], list):
                for row in matrix:
                    degree = tuple(row[:-1])
                    mult = row[-1]
                    results.append((degree, mult))
            elif matrix and isinstance(matrix[0], int):
                # Single row without nesting [a,b,c]
                degree = tuple(matrix[:-1])
                mult = matrix[-1]
                results.append((degree, mult))
    return results


def classify_operators(pe_raw, rcharges_per_type, rep_type_names):
    """Classify GIOs by R-charge.

    Args:
        pe_raw: list of (degree_tuple, multiplicity) from PE
        rcharges_per_type: list of R-charges, one per rep type (gauge-charged only)
        rep_type_names: list of rep type names for display

    Returns:
        dict with categorized operators
    """
    all_ops = []
    relevant = []      # R < 2
    super_relevant = []  # R < 4/3
    unitarity_viol = []  # R < 2/3

    for degree, mult in pe_raw:
        if all(d == 0 for d in degree):
            continue  # skip trivial (identity)

        # R-charge = sum(n_i * r_i)
        R = sum(d * r for d, r in zip(degree, rcharges_per_type))

        op_info = {
            'degree': list(degree),
            'multiplicity': mult,
            'R': float(R),
            'description': ' '.join(
                f"{rep_type_names[i]}^{degree[i]}"
                for i in range(len(degree)) if degree[i] > 0
            ),
        }
        all_ops.append(op_info)

        if R < 2/3:
            unitarity_viol.append(op_info)
        if R < 4/3:
            super_relevant.append(op_info)
        if R < 2:
            relevant.append(op_info)

    return {
        'all_operators': sorted(all_ops, key=lambda x: x['R']),
        'relevant': sorted(relevant, key=lambda x: x['R']),
        'super_relevant': sorted(super_relevant, key=lambda x: x['R']),
        'unitarity_violating': sorted(unitarity_viol, key=lambda x: x['R']),
    }


def analyze_seed(group, seed_entry, all_reps, order=6):
    """Full analysis of a single seed theory.

    Args:
        group: LiE group name (e.g. "A1", "G2")
        seed_entry: entry from consistent_seeds.json
        all_reps: data from all_reps.json
        order: Hilbert series truncation order

    Returns:
        dict with full analysis
    """
    matter = seed_entry['matter']
    desc = seed_entry['description']
    a_val = seed_entry['a']
    c_val = seed_entry['c']
    rcharges = seed_entry['rcharges']

    # Build rep type list for Hilbert series (gauge-charged only)
    dynkin_labels = []
    multiplicities = []
    rep_type_names = []
    rcharges_per_type = []

    # R-charge naming in filter_seeds.py: r{matter_idx}{copy_idx}
    # e.g. matter = [([1,0], 1, 'real', 7), ([0,1], 1, 'real', 14)]
    # -> fields: r01 (matter_idx=0, copy=1), r11 (matter_idx=1, copy=1)
    # For complex: r01, r0b1 (R and R_bar)

    for matter_idx, (label, n, reality, dim_r) in enumerate(matter):
        # Find R-charge: key pattern is r{matter_idx}1 (copy 1, all copies same for W=0)
        r_key = f"r{matter_idx}1"
        r_charge = None
        for k, v in rcharges.items():
            if k == r_key:
                try:
                    r_charge = float(v.split('`')[0])
                except:
                    pass
                break

        if r_charge is None:
            # Try without copy index
            for k, v in rcharges.items():
                if k.startswith(f"r{matter_idx}") and 'b' not in k:
                    try:
                        r_charge = float(v.split('`')[0])
                    except:
                        pass
                    break

        if r_charge is None:
            print(f"  WARNING: could not find R-charge for matter_idx={matter_idx}, "
                  f"label={label}. Keys: {list(rcharges.keys())}")
            r_charge = 0.5  # fallback

        if reality == 'complex':
            # R type
            dynkin_labels.append(label)
            multiplicities.append(n)
            rep_type_names.append(f"{label}")
            rcharges_per_type.append(r_charge)

            # R_bar type (same R-charge)
            rep_lookup = {tuple(r['dynkin_label']): r for r in all_reps[group]['reps']}
            conj_label = label  # default
            if tuple(label) in rep_lookup and rep_lookup[tuple(label)].get('contragr'):
                conj_label = rep_lookup[tuple(label)]['contragr']
            dynkin_labels.append(conj_label)
            multiplicities.append(n)
            rep_type_names.append(f"{conj_label}b")
            rcharges_per_type.append(r_charge)
        else:
            dynkin_labels.append(label)
            multiplicities.append(n)
            rep_type_names.append(f"{label}")
            rcharges_per_type.append(r_charge)

    # Compute Hilbert series
    hs = compute_hilbert_series(
        group, dynkin_labels, multiplicities, order,
        compute_pl=False
    )

    # Parse and classify
    pe_data = parse_pe_raw(hs['pe_raw'])
    classification = classify_operators(pe_data, rcharges_per_type, rep_type_names)

    return {
        'group': group,
        'description': desc,
        'a': a_val,
        'c': c_val,
        'rcharges': rcharges,
        'rcharges_per_type': {n: r for n, r in zip(rep_type_names, rcharges_per_type)},
        'n_relevant': len(classification['relevant']),
        'n_super_relevant': len(classification['super_relevant']),
        'n_unitarity_violating': len(classification['unitarity_violating']),
        'relevant_ops': classification['relevant'],
        'super_relevant_ops': classification['super_relevant'],
        'unitarity_violating_ops': classification['unitarity_violating'],
        'hs_time': hs['time_total'],
    }


if __name__ == '__main__':
    with open('data/all_reps.json') as f:
        all_reps = json.load(f)
    with open('data/consistent_seeds.json') as f:
        seeds = json.load(f)

    groups = sys.argv[1:] if len(sys.argv) > 1 else ['G2']
    order = 6

    for group in groups:
        if group not in seeds:
            print(f"No seeds for {group}")
            continue

        print(f"\n{'='*70}")
        print(f"  {group}: {len(seeds[group])} consistent seeds")
        print(f"{'='*70}")

        for i, seed in enumerate(seeds[group][:5]):  # first 5 for testing
            print(f"\n  --- Seed {i+1}: {seed['description']} ---")
            result = analyze_seed(group, seed, all_reps, order=order)
            print(f"  a={result['a'][:8]}, c={result['c'][:8]}")
            print(f"  Relevant: {result['n_relevant']}, "
                  f"Super-relevant: {result['n_super_relevant']}, "
                  f"Unitarity-violating: {result['n_unitarity_violating']}")
            if result['relevant_ops']:
                print(f"  Relevant operators:")
                for op in result['relevant_ops'][:5]:
                    print(f"    R={op['R']:.4f}  mult={op['multiplicity']}  {op['description']}")
