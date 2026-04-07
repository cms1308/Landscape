"""
Enumerate all asymptotically free matter content for each gauge group.

Conditions:
1. Asymptotic freedom: b_0 = 3*h^v - sum_i N_i * T(R_i) > 0
2. Anomaly cancellation: sum_i N_i * A(R_i) = 0
   - For real/pseudo-real reps: A(R) = 0 automatically
   - For complex reps: A(R) = -A(R_bar), so need to balance R and R_bar
   - Only SU(N) for N >= 3 and E6 have nonzero cubic anomaly
   - For our groups (A1-A4, B, C, D, G2, F4): only A2, A3, A4 need anomaly cancellation

Strategy:
- For groups with no cubic anomaly (A1, B, C, D, G2, F4):
  enumerate all multisets of reps with total T < 3*h^v
- For A_r (r >= 2):
  complex reps must appear with conjugates; the cubic anomaly coefficient
  A(R) for SU(N) is computed from the Dynkin labels.
  For simplicity, we require complex reps to come in (R, R_bar) pairs.
  This automatically cancels the cubic anomaly since A(R_bar) = -A(R).
"""

import json
from fractions import Fraction
from itertools import product as iprod


def load_reps(filename='data/all_reps.json'):
    with open(filename) as f:
        return json.load(f)


def has_cubic_anomaly(group):
    """Groups with potentially nonzero d^{abc} symbol."""
    # SU(N) for N >= 3 (i.e. A_r for r >= 2) and E6
    return group.startswith('A') and int(group[1:]) >= 2 or group == 'E6'


def enumerate_matter(group_data):
    """Enumerate all matter content with b_0 > 0 and anomaly cancellation.

    Returns list of matter content, each as a list of (dynkin_label, N_copies).
    For complex reps, the conjugate is implicitly included.
    """
    group = group_data['group']
    h_v = group_data['h_dual']
    T_budget = Fraction(3 * h_v)  # b_0 > 0 means total T < 3*h^v
    reps = group_data['reps']

    # Build list of "slots" to fill:
    # For real/pseudo-real reps: each copy contributes T(R) to the budget
    # For complex reps: each copy means (R + R_bar), contributing 2*T(R) to the budget
    #   (since T(R) = T(R_bar) and we must include both for anomaly cancellation)
    slots = []
    for r in reps:
        T_R = Fraction(r['T_R'])
        is_complex = r['reality'] == 'complex'
        label = tuple(r['dynkin_label'])
        conj = tuple(r['contragr']) if r['contragr'] else None

        if is_complex:
            # Each "copy" means one R + one R_bar, costing 2*T(R)
            slots.append({
                'label': label,
                'conj': conj,
                'T_per_copy': 2 * T_R,  # R + R_bar
                'reality': 'complex',
                'dim': r['dim'],
            })
        else:
            slots.append({
                'label': label,
                'conj': None,
                'T_per_copy': T_R,
                'reality': r['reality'],
                'dim': r['dim'],
            })

    # Enumerate all multisets: for each slot, choose N_i = 0, 1, 2, ...
    # such that sum N_i * T_per_copy < T_budget
    # Use recursive enumeration with pruning
    results = []

    def search(idx, remaining_T, current):
        """Recursively enumerate valid matter content."""
        # Record current combination if non-empty
        if current:
            results.append(list(current))

        # Try adding more copies of slots[idx:]
        for i in range(idx, len(slots)):
            T_cost = slots[i]['T_per_copy']
            if T_cost <= 0:
                continue
            max_copies = int(remaining_T / T_cost)  # floor
            if remaining_T / T_cost == max_copies and max_copies > 0:
                max_copies -= 1  # strict inequality b_0 > 0

            for n in range(1, max_copies + 1):
                new_T = remaining_T - n * T_cost
                if new_T <= 0:
                    break
                current.append((slots[i]['label'], n, slots[i]))
                search(i + 1, new_T, current)
                current.pop()

    search(0, T_budget, [])
    return results


def format_matter(matter_list):
    """Format matter content for display."""
    parts = []
    for label, n, slot in matter_list:
        if slot['reality'] == 'complex':
            if n == 1:
                parts.append(f"{list(label)}+conj")
            else:
                parts.append(f"{n}x({list(label)}+conj)")
        else:
            if n == 1:
                parts.append(f"{list(label)}")
            else:
                parts.append(f"{n}x{list(label)}")
    return ' + '.join(parts)


def compute_b0(matter_list, h_v):
    """Compute b_0 for a matter content."""
    total_T = sum(Fraction(n) * slot['T_per_copy'] for _, n, slot in matter_list)
    return Fraction(3 * h_v) - total_T


if __name__ == '__main__':
    import sys
    all_reps = load_reps()

    groups = sys.argv[1:] if len(sys.argv) > 1 else sorted(all_reps.keys())

    all_seeds = {}
    for grp in groups:
        if grp not in all_reps:
            continue
        gdata = all_reps[grp]

        matters = enumerate_matter(gdata)

        print(f"\n{'='*80}")
        print(f"  {grp}  |  h^v = {gdata['h_dual']}  |  3*h^v = {3*gdata['h_dual']}  |  {len(matters)} matter contents")
        print(f"{'='*80}")

        for m in sorted(matters, key=lambda x: sum(Fraction(n) * s['T_per_copy'] for _, n, s in x)):
            b0 = compute_b0(m, gdata['h_dual'])
            desc = format_matter(m)
            total_T = Fraction(3 * gdata['h_dual']) - b0
            print(f"  b0={str(b0):>6}  T_total={str(total_T):>6}  {desc}")

        all_seeds[grp] = {
            'h_dual': gdata['h_dual'],
            'T_budget': 3 * gdata['h_dual'],
            'n_seeds': len(matters),
            'seeds': [
                {
                    'matter': [(list(label), n, slot['reality'], slot['dim'])
                               for label, n, slot in m],
                    'b0': str(compute_b0(m, gdata['h_dual'])),
                    'description': format_matter(m),
                }
                for m in matters
            ],
        }

    with open('data/seed_matters.json', 'w') as f:
        json.dump(all_seeds, f, indent=2)
    print(f"\nSaved to data/seed_matters.json")
