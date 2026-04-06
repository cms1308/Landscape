"""
Compute representation data for simple Lie groups using LiE.

The Dynkin index T(R) is computed via:
  T(R) = dim(R) * C_2(R) / dim(G)

where C_2(R) = (lambda, lambda + 2*rho) with the inner product normalized
so that (alpha_long, alpha_long) = 2.

In Dynkin label coordinates, the inner product is:
  (lambda, mu) = sum_{ij} lambda_i * D_i * (A^{-1})_{ij} * mu_j

where A is the Cartan matrix and D = diag(d_1,...,d_r) with d_i = 2/||alpha_i||^2.
For simply-laced (A,D,E): d_i = 1 for all i.
For B_r: d_i = 1 for i < r, d_r = 2.
For C_r: d_i = 1 for i < r, d_r = 1/2.  (Wait, C has short roots at position r?)
Actually: for C_r, the LONG root is alpha_1,...,alpha_{r-1} with ||alpha||^2 = 2,
and alpha_r is long with ||alpha_r||^2 = 4... no.

Standard Bourbaki: for C_r, alpha_r is the long root.
d_i = 2/||alpha_i||^2:
  B_r: alpha_1,...,alpha_{r-1} long (||^2=2, d=1), alpha_r short (||^2=1, d=2)
  C_r: alpha_1,...,alpha_{r-1} short (||^2=1, d=2), alpha_r long (||^2=2, d=1)
  G_2: alpha_1 short (d=3), alpha_2 long (d=1)
  F_4: alpha_1,alpha_2 long (d=1), alpha_3,alpha_4 short (d=2)

Actually LiE uses its own conventions. Let me just use LiE to compute everything
self-consistently via tensor product decomposition.

Method: T(R) can be extracted from R ⊗ R* = adj + ... + trivial.
  Tr_R(T^a T^b) = T(R) delta^{ab}
  => sum_R (T^a)_{ij} (T^b)_{ji} = T(R) delta^{ab}
  => Tr_{R⊗R*}(T^a ⊗ 1 + 1 ⊗ T^a)(T^b ⊗ 1 + 1 ⊗ T^b) includes cross terms...

Actually simplest: use the identity
  dim(R) * C_2(R) = T(R) * dim(G)
and compute C_2(R) properly.

FINAL APPROACH: compute everything in Python using the inverse Cartan matrix
and the symmetrization factors d_i, obtained from LiE's Cartan matrix.
"""

import subprocess
import json
from fractions import Fraction

# Symmetrization factors d_i = diagonal of D where D*A = symmetrized Cartan matrix
# d_i can be read from: d_i * A_{ii} = 2 (diagonal of symmetrized matrix)
# => d_i = 2 / A_{ii}  ... NO, that's only if off-diag is symmetric
# Actually: D*A is symmetric. D = diag(d_1,...,d_r).
# (D*A)_{ij} = d_i * A_{ij} should equal (D*A)_{ji} = d_j * A_{ji}
# So d_i * A_{ij} = d_j * A_{ji}
# Since A_{ii} = 2: d_i * 2 = d_i * A_{ii}, that's trivially true.
# From A_{ij}/A_{ji} = d_j/d_i when i≠j, A_{ij}≠0.

# For our purposes, we use the KNOWN d_i values:
D_FACTORS = {
    'A': None,  # all 1
    'B': None,  # computed: 1,...,1,2
    'C': None,  # computed: 2,...,2,1 (LiE convention)
    'D': None,  # all 1
    'E': None,  # all 1
    'F': None,  # 1,1,2,2 (LiE convention for F4)
    'G': None,  # 3,1 (LiE convention for G2)
}


def lie_eval(code):
    """Run LiE code and return stdout lines."""
    res = subprocess.run('lie', input=code, capture_output=True, encoding='UTF-8', shell=True)
    lines = [l.strip() for l in res.stdout.split('\n') if l.strip()]
    return lines


def get_d_factors(group):
    """Get symmetrization factors d_i from Cartan matrix."""
    code = f"setdefault {group}\nprint(Cartan())\nprint(Lie_rank)\n"
    lines = lie_eval(code)
    # Parse rank
    rank = int(lines[-1])
    # Parse Cartan matrix
    cartan_str = ' '.join(lines[:-1])
    # Extract numbers
    import re
    nums = [int(x) for x in re.findall(r'-?\d+', cartan_str)]
    A = []
    for i in range(rank):
        A.append(nums[i*rank:(i+1)*rank])

    # Compute d_i = 2/||alpha_i||^2 using the relation d_i * A_{ij} = d_j * A_{ji}
    # Start with d_0 = 1 and propagate (ratios are correct)
    d = [Fraction(1)] * rank
    for i in range(rank):
        for j in range(i+1, rank):
            if A[i][j] != 0:
                d[j] = d[i] * Fraction(A[i][j], A[j][i])
                break

    # Normalize: the long root has ||alpha||^2 = 2, so d_long = 1.
    # The long root has the SMALLEST d_i (since d = 2/||alpha||^2 and long root has largest ||alpha||^2).
    # So normalize by dividing all d by min(d).
    d_min = min(d)
    d = [x / d_min for x in d]

    return rank, A, d


def compute_T(rank, A, d, dim_G, lam):
    """Compute Dynkin index T(R) for representation with Dynkin label lam.

    C_2(R) = sum_{ij} lam_i * d_i * (A^{-1})_{ij} * (lam_j + 2)
    T(R) = dim(R) * C_2(R) / dim(G)

    We compute A^{-1} via the adjugate: A^{-1} = adj(A) / det(A)
    """
    # Compute A inverse using fractions for exact arithmetic
    from fractions import Fraction
    # Convert A to Fraction matrix
    AF = [[Fraction(A[i][j]) for j in range(rank)] for i in range(rank)]

    # Gaussian elimination to find inverse
    # Augment [A | I]
    aug = [AF[i] + [Fraction(1) if i == j else Fraction(0) for j in range(rank)] for i in range(rank)]
    for col in range(rank):
        # Find pivot
        pivot = None
        for row in range(col, rank):
            if aug[row][col] != 0:
                pivot = row
                break
        aug[col], aug[pivot] = aug[pivot], aug[col]
        # Scale
        scale = aug[col][col]
        aug[col] = [x / scale for x in aug[col]]
        # Eliminate
        for row in range(rank):
            if row != col and aug[row][col] != 0:
                factor = aug[row][col]
                aug[row] = [aug[row][j] - factor * aug[col][j] for j in range(2*rank)]
    Ainv = [row[rank:] for row in aug]

    # rho in Dynkin labels = [1,1,...,1]
    lam_plus_2rho = [lam[i] + 2 for i in range(rank)]

    # Metric g = A^{-1} . D where D_j = ||alpha_j||^2 / 2 = 1/d_j
    # (d_j = 2/||alpha_j||^2 as computed from Cartan matrix)
    # g_{ij} = sum_k Ainv_{ik} * D_k * delta_{kj} = Ainv_{ij} * (1/d_j)
    # C_2 = sum_{ij} lam_i * g_{ij} * (lam_j + 2)
    #      = sum_{ij} lam_i * Ainv_{ij} * (1/d_j) * (lam_j + 2)
    C2 = Fraction(0)
    for i in range(rank):
        for j in range(rank):
            C2 += Fraction(lam[i]) * Ainv[i][j] / d[j] * Fraction(lam_plus_2rho[j])

    # Get dimension from LiE
    # (passed in from caller)
    return C2


def get_dim(group, label_str):
    code = f"setdefault {group}\nprint(dim({label_str}))\n"
    lines = lie_eval(code)
    return int(lines[0])


def get_reality(group, label_str):
    """Return 'real', 'pseudo-real', or 'complex'."""
    code = f"setdefault {group}\nx = {label_str}\nprint(contragr(x) == x)\n"
    lines = lie_eval(code)
    is_self_conj = int(lines[0])
    if not is_self_conj:
        return 'complex'

    code2 = f"setdefault {group}\nx = {label_str}\nz = null(Lie_rank)\na = alt_tensor(2,x)\nprint(a|z)\n"
    lines2 = lie_eval(code2)
    trivial_in_alt2 = int(lines2[0])
    if trivial_in_alt2 > 0:
        return 'pseudo-real'
    else:
        return 'real'


def get_adj_label(group):
    code = f"setdefault {group}\nprint(adjoint())\n"
    lines = lie_eval(code)
    s = lines[0]
    bracket = s[s.index('['):s.index(']')+1]
    import ast
    return ast.literal_eval(bracket)


# Groups and their representations
GROUPS_REPS = {
    'A1': [
        ('fund', [1]),
        ('adj', [2]),
        ('[3]', [3]),
        ('[4]', [4]),
    ],
    'A2': [
        ('fund', [1,0]),
        ('anti-fund', [0,1]),
        ('adj', [1,1]),
        ('sym', [2,0]),
        ('sym-bar', [0,2]),
        ('[3,0]', [3,0]),
    ],
    'A3': [
        ('fund', [1,0,0]),
        ('anti-fund', [0,0,1]),
        ('anti-sym', [0,1,0]),
        ('adj', [1,0,1]),
        ('sym', [2,0,0]),
        ('sym-bar', [0,0,2]),
    ],
    'A4': [
        ('fund', [1,0,0,0]),
        ('anti-fund', [0,0,0,1]),
        ('anti-sym', [0,1,0,0]),
        ('anti-sym-bar', [0,0,1,0]),
        ('adj', [1,0,0,1]),
        ('sym', [2,0,0,0]),
    ],
    'B2': [
        ('vector', [1,0]),
        ('spinor', [0,1]),
        ('adj', [0,2]),
        ('sym traceless', [2,0]),
    ],
    'C2': [
        ('fund', [1,0]),
        ('[0,1]', [0,1]),
        ('adj', [2,0]),
        ('[0,2]', [0,2]),
        ('[1,1]', [1,1]),
    ],
    'G2': [
        ('fund', [1,0]),
        ('adj', [0,1]),
        ('[2,0]', [2,0]),
        ('[1,1]', [1,1]),
    ],
    'B3': [
        ('vector', [1,0,0]),
        ('spinor', [0,0,1]),
        ('adj', [0,1,0]),
        ('sym traceless', [2,0,0]),
    ],
    'C3': [
        ('fund', [1,0,0]),
        ('[0,1,0]', [0,1,0]),
        ('adj', [2,0,0]),
        ('[0,0,1]', [0,0,1]),
    ],
    'D4': [
        ('vector', [1,0,0,0]),
        ('spinor', [0,0,1,0]),
        ('conj-spinor', [0,0,0,1]),
        ('adj', [0,1,0,0]),
    ],
    'F4': [
        ('fund', [0,0,0,1]),
        ('adj', [1,0,0,0]),
    ],
}


if __name__ == '__main__':
    all_results = {}

    for grp, reps in GROUPS_REPS.items():
        rank, A, d = get_d_factors(grp)
        adj_label = get_adj_label(grp)
        dim_G = get_dim(grp, str(adj_label))

        # Verify T(adj) = h^v by computing C2(adj)
        C2_adj = compute_T(rank, A, d, dim_G, adj_label)
        T_adj = C2_adj * dim_G / (2 * dim_G)  # T(adj) = C2(adj) / 2 in physics convention

        print(f"\n{'='*75}")
        print(f"  {grp}  |  dim(G) = {dim_G}  |  d = {[str(x) for x in d]}  |  T(adj) = {T_adj}")
        print(f"{'='*75}")
        print(f"  {'Name':<20} {'Dynkin':<18} {'dim':>5} {'T(R)':>12} {'C2(R)':>12} {'Reality':<12}")
        print(f"  {'-'*75}")

        grp_results = []
        for name, label in reps:
            label_str = str(label)
            dim_R = get_dim(grp, label_str)
            C2_R = compute_T(rank, A, d, dim_G, label)
            T_R = C2_R * dim_R / (2 * dim_G)
            reality = get_reality(grp, label_str)

            print(f"  {name:<20} {label_str:<18} {dim_R:>5} {str(T_R):>12} {str(C2_R):>12} {reality:<12}")
            grp_results.append({
                'name': name,
                'dynkin_label': label,
                'dim': dim_R,
                'T_R': str(T_R),
                'C2_R': str(C2_R),
                'reality': reality,
            })

        all_results[grp] = {
            'dim_G': dim_G,
            'T_adj': str(T_adj),
            'd_factors': [str(x) for x in d],
            'reps': grp_results,
        }

    with open('data/rep_data.json', 'w') as f:
        json.dump(all_results, f, indent=2)
    print(f"\nSaved to data/rep_data.json")
