"""
Systematically enumerate all representations R with T(R) < 3*h^v for each simple Lie group.
This is the single-copy asymptotic freedom bound: a theory with one copy of R
(and nothing else) has b_0 = 3*T(adj) - T(R) = 3*h^v - T(R) > 0.

Strategy: iterate over Dynkin labels [l1, l2, ..., lr] with increasing values,
compute T(R) via the weight space inner product, and keep those with T(R) < 3*h^v.

We bound the search by noting that T(R) grows with the Dynkin labels, so we can
prune branches where T is already too large.
"""

import subprocess
import json
from fractions import Fraction
from itertools import product

# Dual Coxeter numbers
DUAL_COXETER = {
    'A1': 2, 'A2': 3, 'A3': 4, 'A4': 5, 'A5': 6,
    'B2': 3, 'B3': 5, 'B4': 7,
    'C2': 3, 'C3': 4, 'C4': 5,
    'D4': 6, 'D5': 8,
    'G2': 4, 'F4': 9, 'E6': 12, 'E7': 18, 'E8': 30,
}


def lie_eval(code):
    res = subprocess.run('lie', input=code, capture_output=True, encoding='UTF-8', shell=True)
    return [l.strip() for l in res.stdout.split('\n') if l.strip()]


def get_cartan_and_rank(group):
    import re
    code = f"setdefault {group}\nprint(Cartan())\nprint(Lie_rank)\n"
    lines = lie_eval(code)
    rank = int(lines[-1])
    cartan_str = ' '.join(lines[:-1])
    nums = [int(x) for x in re.findall(r'-?\d+', cartan_str)]
    A = []
    for i in range(rank):
        A.append(nums[i*rank:(i+1)*rank])
    return rank, A


def compute_d_factors(rank, A):
    d = [Fraction(1)] * rank
    for i in range(rank):
        for j in range(i+1, rank):
            if A[i][j] != 0:
                d[j] = d[i] * Fraction(A[i][j], A[j][i])
                break
    d_min = min(d)
    d = [x / d_min for x in d]
    return d


def compute_Ainv(rank, A):
    AF = [[Fraction(A[i][j]) for j in range(rank)] for i in range(rank)]
    aug = [AF[i] + [Fraction(1) if i == j else Fraction(0) for j in range(rank)] for i in range(rank)]
    for col in range(rank):
        pivot = None
        for row in range(col, rank):
            if aug[row][col] != 0:
                pivot = row
                break
        aug[col], aug[pivot] = aug[pivot], aug[col]
        scale = aug[col][col]
        aug[col] = [x / scale for x in aug[col]]
        for row in range(rank):
            if row != col and aug[row][col] != 0:
                factor = aug[row][col]
                aug[row] = [aug[row][j] - factor * aug[col][j] for j in range(2*rank)]
    return [row[rank:] for row in aug]


def compute_T(rank, Ainv, d, dim_G, lam):
    """Compute T(R) = dim(R) * C2(R) / (2*dim(G)) where C2 uses metric g = A^{-1}.D"""
    lam_plus_2rho = [lam[i] + 2 for i in range(rank)]
    C2 = Fraction(0)
    for i in range(rank):
        for j in range(rank):
            C2 += Fraction(lam[i]) * Ainv[i][j] / d[j] * Fraction(lam_plus_2rho[j])
    return C2


def get_dim_lie(group, label):
    label_str = '[' + ','.join(str(x) for x in label) + ']'
    code = f"setdefault {group}\nprint(dim({label_str}))\n"
    lines = lie_eval(code)
    return int(lines[0])


def get_contragr_lie(group, label):
    label_str = '[' + ','.join(str(x) for x in label) + ']'
    code = f"setdefault {group}\nprint(contragr({label_str}))\n"
    lines = lie_eval(code)
    import re
    nums = [int(x) for x in re.findall(r'-?\d+', lines[0])]
    return tuple(nums)


def get_reality_lie(group, label):
    label_str = '[' + ','.join(str(x) for x in label) + ']'
    code = f"setdefault {group}\nx = {label_str}\nprint(contragr(x) == x)\n"
    lines = lie_eval(code)
    is_self_conj = int(lines[0])
    if not is_self_conj:
        return 'complex'
    code2 = f"setdefault {group}\nx = {label_str}\nz = null(Lie_rank)\na = alt_tensor(2,x)\nprint(a|z)\n"
    lines2 = lie_eval(code2)
    trivial_in_alt2 = int(lines2[0])
    return 'pseudo-real' if trivial_in_alt2 > 0 else 'real'


def get_adj_label(group):
    code = f"setdefault {group}\nprint(adjoint())\n"
    lines = lie_eval(code)
    s = lines[0]
    bracket = s[s.index('['):s.index(']')+1]
    import ast
    return list(ast.literal_eval(bracket))


def enumerate_reps(group):
    """Enumerate all Dynkin labels with T(R) < 3*h^v."""
    h_v = DUAL_COXETER[group]
    T_max = Fraction(3 * h_v)

    rank, A = get_cartan_and_rank(group)
    d = compute_d_factors(rank, A)
    Ainv = compute_Ainv(rank, A)
    adj_label = get_adj_label(group)
    dim_G = get_dim_lie(group, adj_label)

    # We need to search over Dynkin labels [l1, ..., lr] with l_i >= 0.
    # T(R) grows roughly quadratically with the labels, so we can bound the search.
    # For a single label l_i with others=0: C2 ~ l_i*(l_i+2) * Ainv[i][i] / d[i]
    # and dim grows, so T grows fast. We search up to some max per component.
    # Use max_label = 20 as a safe upper bound (will prune via T).
    max_label = 20

    results = []
    seen_contragr = set()  # track conjugate pairs to avoid duplicates

    # Generate all labels with bounded components
    # For efficiency, use recursive generation with pruning
    def search(pos, label_so_far):
        if pos == rank:
            lam = list(label_so_far)
            if all(l == 0 for l in lam):
                return  # skip trivial rep

            # Compute C2 (don't need dim yet for the bound check)
            C2 = compute_T(rank, Ainv, d, dim_G, lam)

            # Quick bound: C2 must be positive for non-trivial reps
            if C2 <= 0:
                return

            # Get dimension
            dim_R = get_dim_lie(group, lam)
            T_R = C2 * dim_R / (2 * dim_G)

            if T_R < T_max:
                lam_tuple = tuple(lam)
                # Check if contragredient already seen
                contragr = get_contragr_lie(group, lam)
                if contragr < lam_tuple and contragr in seen_contragr:
                    # Already have the conjugate, skip (will be listed as conjugate)
                    pass
                else:
                    reality = get_reality_lie(group, lam)
                    results.append({
                        'dynkin_label': lam,
                        'dim': dim_R,
                        'T_R': str(T_R),
                        'T_R_float': float(T_R),
                        'C2': str(C2),
                        'reality': reality,
                        'contragr': list(contragr) if contragr != lam_tuple else None,
                    })
                    seen_contragr.add(lam_tuple)
            return

        # Try values for position pos
        for val in range(0, max_label + 1):
            new_label = label_so_far + (val,)
            # Partial C2 lower bound: even with remaining labels = 0,
            # C2 is at least what we have so far. But computing partial C2
            # is complex, so just use a simple bound on the label sum.
            # For very high labels, T will be huge, so we can prune.
            # Simple heuristic: if any single label > 2*max_label_estimate, skip.
            if val > max_label:
                break
            search(pos + 1, new_label)

    # For efficiency with high-rank groups, limit the search space
    # Estimate: for each position, the max useful label is roughly
    # sqrt(3*h^v * 2*dim_G / (Ainv[i][i]/d[i])) but this is complex.
    # Just use a generous bound and let T pruning handle it.
    # For rank > 4, reduce max_label to avoid combinatorial explosion.
    if rank <= 2:
        max_label = 20
    elif rank <= 4:
        max_label = 10
    elif rank <= 6:
        max_label = 6
    else:
        max_label = 4

    # Instead of recursive search (too slow for high rank), use iterative approach
    # with early termination based on C2 partial sums
    print(f"\n{'='*80}")
    print(f"  {group}  |  rank = {rank}  |  dim(G) = {dim_G}  |  h^v = {h_v}  |  T_max = {3*h_v}")
    print(f"  d = {[str(x) for x in d]}")
    print(f"{'='*80}")

    # Generate candidates: use itertools.product with bounded range
    # For rank 1-2: exhaustive. For higher rank: be smarter.
    from itertools import product as iprod

    total_checked = 0
    for lam_tuple in iprod(*(range(0, max_label+1) for _ in range(rank))):
        lam = list(lam_tuple)
        if all(l == 0 for l in lam):
            continue

        # Compute C2
        lam_plus_2rho = [lam[i] + 2 for i in range(rank)]
        C2 = Fraction(0)
        for i in range(rank):
            for j in range(rank):
                C2 += Fraction(lam[i]) * Ainv[i][j] / d[j] * Fraction(lam_plus_2rho[j])

        if C2 <= 0:
            continue

        # Quick upper bound on T: T = C2 * dim / (2*dim_G)
        # dim >= 1, so if C2 / (2*dim_G) >= T_max already, skip (but dim >> 1 so this rarely helps)
        # Actually dim grows fast, so we need the actual dim.
        # But calling LiE for every candidate is slow. Let's compute dim only for promising C2.
        # Heuristic: if C2 > T_max * 2 * dim_G (which means T > T_max even for dim=1), skip.
        # But dim >= max(l_i)+1 roughly, so a tighter bound: skip if C2 > T_max * 2 * dim_G.
        if C2 > T_max * 2 * dim_G:
            continue

        dim_R = get_dim_lie(group, lam)
        T_R = C2 * dim_R / (2 * dim_G)
        total_checked += 1

        if T_R < T_max:
            contragr = get_contragr_lie(group, lam)
            # For complex reps, keep the one with LOWER Dynkin label tuple
            # (e.g. [1,0] before [0,1] for A2) and list the other as conjugate
            if tuple(contragr) != lam_tuple:
                # Complex rep — check if conjugate already seen or if we should skip
                if tuple(contragr) in seen_contragr:
                    continue  # conjugate already listed as primary
                if lam_tuple < tuple(contragr):
                    continue  # we're the lower label; wait for the higher one
            reality = get_reality_lie(group, lam)
            results.append({
                'dynkin_label': lam,
                'dim': dim_R,
                'T_R': str(T_R),
                'T_R_float': float(T_R),
                'C2': str(C2),
                'reality': reality,
                'contragr': list(contragr) if tuple(contragr) != lam_tuple else None,
            })
            seen_contragr.add(lam_tuple)

    # Sort by T(R)
    results_sorted = sorted(results, key=lambda x: x['T_R_float'])

    print(f"  Checked {total_checked} candidates, found {len(results_sorted)} reps with T < {3*h_v}")
    print(f"\n  {'Dynkin':<20} {'dim':>5} {'T(R)':>12} {'Reality':<12} {'Conjugate':<20}")
    print(f"  {'-'*70}")
    for r in results_sorted:
        conj_str = str(r['contragr']) if r['contragr'] else 'self'
        print(f"  {str(r['dynkin_label']):<20} {r['dim']:>5} {r['T_R']:>12} {r['reality']:<12} {conj_str:<20}")

    return {
        'group': group,
        'rank': rank,
        'dim_G': dim_G,
        'h_dual': h_v,
        'T_max': 3 * h_v,
        'd_factors': [str(x) for x in d],
        'reps': results_sorted,
    }


GROUPS = ['A1', 'A2', 'A3', 'A4', 'B2', 'B3', 'C2', 'C3', 'D4', 'G2', 'F4']

if __name__ == '__main__':
    import sys
    groups = sys.argv[1:] if len(sys.argv) > 1 else GROUPS

    all_results = {}
    for grp in groups:
        if grp not in DUAL_COXETER:
            print(f"Unknown group: {grp}")
            continue
        all_results[grp] = enumerate_reps(grp)

    with open('data/all_reps.json', 'w') as f:
        json.dump(all_results, f, indent=2)
    print(f"\nSaved to data/all_reps.json")
