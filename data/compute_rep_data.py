"""
Compute representation data for simple Lie groups using LiE.
For each group and representation (Dynkin label), computes:
  - dimension |R|
  - Dynkin index T(R) (normalized so T(fund) = 1/2 for A_r)
  - reality: real, pseudo-real, or complex

LiE's norm() uses its own normalization. We extract dim(R)*C2_LiE(R) which
is proportional to T(R), then calibrate using the adjoint:
  T(adj) = h^v in our convention
  dim(R)*C2_LiE(R) / dim(adj)*C2_LiE(adj) = T(R) / T(adj)
  => T(R) = h^v * dim(R)*C2_LiE(R) / (dim(adj)*C2_LiE(adj))
"""

import subprocess
import json
import sys

# Known dual Coxeter numbers and adjoint Dynkin labels
DUAL_COXETER = {
    'A1': 2, 'A2': 3, 'A3': 4, 'A4': 5, 'A5': 6, 'A6': 7, 'A7': 8,
    'B2': 3, 'B3': 5, 'B4': 7, 'B5': 9,
    'C2': 3, 'C3': 4, 'C4': 5, 'C5': 6,
    'D3': 4, 'D4': 6, 'D5': 8,
    'G2': 4, 'F4': 9, 'E6': 12, 'E7': 18, 'E8': 30,
}

# Adjoint Dynkin labels (highest weight of adjoint rep)
# A_r: [1,0,...,0,1] for r>=2, [2] for r=1
# B_r: [0,1,0,...,0]
# C_r: [2,0,...,0]
# D_r: [0,1,0,...,0]
# G2: [0,1], F4: [1,0,0,0], E6: [0,0,0,0,0,1] (or [1,0,0,0,0,1])...
# Actually just compute from LiE
def get_adj_label(group):
    """Get adjoint Dynkin label by parsing LiE's adjoint() output."""
    code = f"setdefault {group}\nprint(adjoint())\n"
    lines = lie_eval(code)
    # Output like "1X[1,1]" — extract the bracket part
    s = lines[0]
    bracket = s[s.index('['):s.index(']')+1]
    return bracket


def lie_eval(code):
    """Run LiE code and return stdout lines."""
    res = subprocess.run('lie', input=code, capture_output=True, encoding='UTF-8', shell=True)
    lines = [l.strip() for l in res.stdout.split('\n') if l.strip()]
    return lines


def get_dim_c2_lie(group, label_str):
    """Return (dim, dim*C2_LiE) for a representation."""
    code = f"""setdefault {group}
rho = all_one(Lie_rank)
x = {label_str}
print(dim(x))
print(dim(x) * (norm(x + rho) - norm(rho)))
"""
    lines = lie_eval(code)
    dim_r = int(lines[0])
    dc2 = int(lines[1])
    return dim_r, dc2


def get_reality(group, label_str):
    """Return 'real', 'pseudo-real', or 'complex'."""
    # Step 1: check self-conjugacy
    code1 = f"setdefault {group}\nprint(contragr({label_str}) == {label_str})\n"
    lines1 = lie_eval(code1)
    is_self_conj = int(lines1[0])
    if not is_self_conj:
        return 'complex'

    # Step 2: for self-conjugate reps, check Frobenius-Schur indicator
    # alt^2(R) contains trivial => pseudo-real; sym^2(R) contains trivial => real
    # Dim of trivial component in alt^2: use LiE's restriction notation
    code2 = f"""setdefault {group}
x = {label_str}
a = alt_tensor(2,x)
print(dim(a))
s = sym_tensor(2,x)
d = dim(x)
ad = dim(a)
sd = dim(s)
print(ad)
print(sd)
"""
    lines2 = lie_eval(code2)
    dim_alt2 = int(lines2[1])
    dim_sym2 = int(lines2[2])
    dim_r = int(lines2[0].split()[0]) if lines2 else 0

    # For dim d rep: dim(sym^2) + dim(alt^2) = d*(d+1)/2 + d*(d-1)/2 = d^2
    # If alt^2 contains trivial: d*(d-1)/2 includes a singlet => pseudo-real
    # Easier: for self-conjugate, check if d*(d-1)/2 == dim_alt2 or dim_alt2 includes +1
    # Actually just check: dim(alt^2) vs d(d-1)/2. If dim(alt^2) includes a singlet
    # from the symplectic form, then it's pseudo-real.
    # Simpler: real reps have Frobenius-Schur = +1, pseudo-real have -1.
    # FS indicator = (1/|G|) sum chi(g^2). For Lie groups:
    # FS = +1 if sym^2 contains trivial, -1 if alt^2 contains trivial.
    # But both sym^2 and alt^2 can contain trivial for different reasons...
    # For irreducible self-conjugate: exactly one of sym^2 or alt^2 contains trivial (once).
    # Check which by comparing: dim(sym^2(R)) vs dim(R)*(dim(R)+1)/2
    # If dim(sym^2) = d*(d+1)/2 then no singlet in sym^2 => singlet in alt^2 => pseudo-real
    # If dim(sym^2) = d*(d+1)/2 + extra => NO, that's wrong direction
    # Actually: R tensor R = sym^2(R) + alt^2(R). Singlet in R x R iff self-conjugate (which it is).
    # The singlet lives in either sym^2 or alt^2, indicating real or pseudo-real.
    # Check: if dim_alt2 = d*(d-1)/2, no singlet in alt^2 => singlet in sym^2 => real
    # If dim_alt2 = d*(d-1)/2 + 1... no that's not right either.
    # Let me just count: total dim of R x R = d^2 = dim_sym2 + dim_alt2
    # Number of singlets in R x R (for self-conjugate irrep) = 1
    # This singlet is in sym^2 (real) or alt^2 (pseudo-real)
    # So: if dim_sym2 > d*(d+1)/2 by 1... no, dims of sym^2 and alt^2 as REPRESENTATIONS
    # are not simply d(d+1)/2 etc. because we decompose into irreps of G.

    # Check if trivial rep appears in alt^2(R): use poly|weight notation
    code3 = f"""setdefault {group}
x = {label_str}
z = null(Lie_rank)
a = alt_tensor(2,x)
print(a|z)
"""
    lines3 = lie_eval(code3)
    trivial_in_alt2 = int(lines3[0])
    if trivial_in_alt2 > 0:
        return 'pseudo-real'
    else:
        return 'real'


def compute_group(group, reps):
    """Compute all rep data for a group."""
    h_v = DUAL_COXETER[group]

    # Get adjoint dim*C2 for calibration
    adj_label = get_adj_label(group)
    adj_code = f"""setdefault {group}
rho = all_one(Lie_rank)
print(dim({adj_label}))
print(dim({adj_label}) * (norm({adj_label} + rho) - norm(rho)))
"""
    adj_lines = lie_eval(adj_code)
    dim_adj = int(adj_lines[0])
    dc2_adj = int(adj_lines[1])

    print(f"\n{'='*70}")
    print(f"  {group}  |  dim(G) = {dim_adj}  |  h^v = {h_v}  |  T(adj) = {h_v}")
    print(f"{'='*70}")
    print(f"  {'Name':<20} {'Dynkin':<18} {'dim':>5} {'T(R)':>12} {'Reality':<12}")
    print(f"  {'-'*67}")

    results = []
    for name, label in reps:
        label_str = '[' + ','.join(str(x) for x in label) + ']'
        dim_r, dc2 = get_dim_c2_lie(group, label_str)
        reality = get_reality(group, label_str)

        # T(R) = h^v * dc2 / dc2_adj
        from fractions import Fraction
        T_R = Fraction(h_v * dc2, dc2_adj)

        print(f"  {name:<20} {str(label):<18} {dim_r:>5} {str(T_R):>12} {reality:<12}")
        results.append({
            'name': name,
            'dynkin_label': label,
            'dim': dim_r,
            'T_R': str(T_R),
            'T_R_float': float(T_R),
            'reality': reality,
        })

    return {
        'group': group,
        'dim_G': dim_adj,
        'h_dual': h_v,
        'T_adj': h_v,
        'reps': results,
    }


# Groups and their representations
GROUPS_REPS = {
    'A1': [
        ('fund [1]', [1]),
        ('adj [2]', [2]),
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
        ('anti-sym [0,1,0]', [0,1,0]),
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
        ('adj', [2,0]),
        ('[0,2]', [0,2]),
    ],
    'C2': [
        ('fund', [1,0]),
        ('[0,1]', [0,1]),
        ('adj [2,0]', [2,0]),
        ('[0,2]', [0,2]),
        ('[1,1]', [1,1]),
    ],
    'G2': [
        ('fund [1,0]', [1,0]),
        ('adj [0,1]', [0,1]),
        ('[2,0]', [2,0]),
        ('[1,1]', [1,1]),
        ('[0,2]', [0,2]),
    ],
    'B3': [
        ('vector', [1,0,0]),
        ('spinor', [0,0,1]),
        ('adj [0,1,0]', [0,1,0]),
        ('sym [2,0,0]', [2,0,0]),
    ],
    'C3': [
        ('fund', [1,0,0]),
        ('[0,1,0]', [0,1,0]),
        ('sym [2,0,0]', [2,0,0]),
        ('[0,0,1]', [0,0,1]),
    ],
    'D4': [
        ('vector', [1,0,0,0]),
        ('spinor', [0,0,1,0]),
        ('conj-spinor', [0,0,0,1]),
        ('adj [0,1,0,0]', [0,1,0,0]),
    ],
    'F4': [
        ('[1,0,0,0]', [1,0,0,0]),
        ('[0,0,0,1]', [0,0,0,1]),
        ('[0,1,0,0]', [0,1,0,0]),
    ],
}


if __name__ == '__main__':
    all_results = {}
    for grp, reps in GROUPS_REPS.items():
        all_results[grp] = compute_group(grp, reps)

    with open('data/rep_data.json', 'w') as f:
        json.dump(all_results, f, indent=2, default=str)
    print(f"\nSaved to data/rep_data.json")
