"""
Find primitive gauge-invariant tensor structures via LiE.

For each (group, representation) pair, determines at which degrees
the symmetric/antisymmetric tensor powers contain a gauge singlet.
Also finds mixed-rep invariants (tensor products of different reps).

This is group-independent: LiE handles the representation theory
for any simple Lie group.
"""

import subprocess
import json
from functools import lru_cache


def lie_eval(code):
    """Run LiE code and return stdout lines."""
    res = subprocess.run('lie', input=code, capture_output=True, encoding='UTF-8', shell=True)
    return [l.strip() for l in res.stdout.split('\n') if l.strip()]


@lru_cache(maxsize=256)
def singlet_in_tensor_power(lie_group, dynkin_label_str, degree, symmetry):
    """Check if the trivial rep appears in sym^d(R) or alt^d(R).

    Args:
        lie_group: e.g. "A1", "G2"
        dynkin_label_str: e.g. "[1]", "[1,0]"
        degree: integer >= 2
        symmetry: "sym" or "alt"

    Returns:
        int: multiplicity of trivial rep (0 if none)
    """
    func = "sym_tensor" if symmetry == "sym" else "alt_tensor"
    code = f"""setdefault {lie_group}
x = {dynkin_label_str}
z = null(Lie_rank)
t = {func}({degree}, x)
if t == 0*1X z then print(0) else print(t|z) fi
"""
    lines = lie_eval(code)
    try:
        return int(lines[0])
    except (IndexError, ValueError):
        return 0


@lru_cache(maxsize=256)
def singlet_in_tensor_product(lie_group, label1_str, label2_str):
    """Check if trivial rep appears in R1 tensor R2.

    Returns:
        int: multiplicity of trivial rep
    """
    code = f"""setdefault {lie_group}
z = null(Lie_rank)
t = tensor({label1_str}, {label2_str})
if t == 0*1X z then print(0) else print(t|z) fi
"""
    lines = lie_eval(code)
    try:
        return int(lines[0])
    except (IndexError, ValueError):
        return 0


def find_single_rep_invariants(lie_group, dynkin_label, max_degree=6):
    """Find all degrees at which sym^d(R) or alt^d(R) contains a singlet.

    Args:
        lie_group: LiE group name
        dynkin_label: list of ints
        max_degree: search up to this degree

    Returns:
        list of (degree, symmetry_type, multiplicity) tuples
        symmetry_type: "sym" or "alt"
    """
    label_str = '[' + ','.join(str(x) for x in dynkin_label) + ']'
    results = []

    for d in range(2, max_degree + 1):
        for sym_type in ["alt", "sym"]:
            mult = singlet_in_tensor_power(lie_group, label_str, d, sym_type)
            if mult > 0:
                results.append((d, sym_type, mult))

    return results


def find_mixed_rep_invariants(lie_group, labels):
    """Find which pairs of rep types have a bilinear invariant (R1 x R2 contains trivial).

    Args:
        lie_group: LiE group name
        labels: list of Dynkin label lists

    Returns:
        list of (i, j, multiplicity) where labels[i] x labels[j] has trivial
    """
    results = []
    for i in range(len(labels)):
        for j in range(i, len(labels)):
            l1_str = '[' + ','.join(str(x) for x in labels[i]) + ']'
            l2_str = '[' + ','.join(str(x) for x in labels[j]) + ']'
            mult = singlet_in_tensor_product(lie_group, l1_str, l2_str)
            if mult > 0:
                results.append((i, j, mult))
    return results


def get_contraction_type(lie_group, dynkin_label):
    """Determine if the bilinear invariant for R x R is symmetric or antisymmetric.

    Returns:
        "sym" if singlet in Sym^2(R), "alt" if singlet in Alt^2(R),
        "both" if both, "none" if neither (R x R has no singlet, i.e. complex rep)
    """
    label_str = '[' + ','.join(str(x) for x in dynkin_label) + ']'
    in_sym = singlet_in_tensor_power(lie_group, label_str, 2, "sym")
    in_alt = singlet_in_tensor_power(lie_group, label_str, 2, "alt")

    if in_sym > 0 and in_alt > 0:
        return "both"
    elif in_sym > 0:
        return "sym"
    elif in_alt > 0:
        return "alt"
    else:
        return "none"


def build_invariant_catalog(lie_group, dynkin_labels, max_degree=6):
    """Build complete catalog of gauge-invariant structures for a set of rep types.

    Args:
        lie_group: LiE group name
        dynkin_labels: list of Dynkin label lists (one per rep type, excluding singlets)
        max_degree: max tensor power to check

    Returns:
        dict with:
          'single_rep': {i: [(degree, sym_type, mult), ...]}
          'mixed_rep': [(i, j, mult), ...]
          'contraction_types': {i: "sym"|"alt"|"both"|"none"}
    """
    single = {}
    for i, label in enumerate(dynkin_labels):
        invs = find_single_rep_invariants(lie_group, label, max_degree)
        if invs:
            single[i] = invs

    mixed = find_mixed_rep_invariants(lie_group, dynkin_labels)

    contraction = {}
    for i, label in enumerate(dynkin_labels):
        contraction[i] = get_contraction_type(lie_group, label)

    return {
        'single_rep': single,
        'mixed_rep': mixed,
        'contraction_types': contraction,
    }


@lru_cache(maxsize=1024)
def singlet_in_multi_tensor(lie_group, *label_strs):
    """Check if trivial rep appears in the tensor product of multiple reps.

    Uses LiE's tensor(vec, vec, grp) iteratively: decompose R1⊗R2 into
    components, then tensor each component with R3, etc. Count total singlets.

    Args:
        lie_group: e.g. "A3"
        *label_strs: Dynkin labels as strings, e.g. "[1,0,0]", "[0,0,1]", "[1,0,1]"

    Returns:
        int: multiplicity of trivial rep in the tensor product
    """
    if len(label_strs) == 0:
        return 1
    if len(label_strs) == 1:
        code = f"setdefault {lie_group}\nz = null(Lie_rank)\nprint(1X{label_strs[0]}|z)\n"
        lines = lie_eval(code)
        try:
            return int(lines[0])
        except (IndexError, ValueError):
            return 0

    # Iterative tensor product via LiE.
    # LiE's tensor(vec, vec, grp) decomposes R1⊗R2 into irreps.
    # For 3+ reps: decompose R1⊗R2, extract each irrep component,
    # tensor each with R3, sum up, repeat.
    # Use LiE's dom_chars to extract components from a polynomial.

    # Build LiE code that iterates:
    # Step 1: t = tensor(R1, R2, G)  → polynomial
    # Step 2: for each irrep V in t with coefficient c:
    #           accumulate c * tensor(V, R3, G)
    # Step 3: repeat with R4, etc.
    # Final: count singlets via |z

    code = f"setdefault {lie_group}\nz = null(Lie_rank)\n"

    if len(label_strs) == 2:
        code += f"t = tensor({label_strs[0]}, {label_strs[1]}, {lie_group})\n"
        code += "print(t|z)\n"
    else:
        # Start with tensor of first two
        code += f"t = tensor({label_strs[0]}, {label_strs[1]}, {lie_group})\n"
        # For each additional rep, expand using expon(poly, j) and length(poly)
        for i in range(2, len(label_strs)):
            code += f"""new_t = 0*1X z
for j = 1 to length(t) do
  new_t = new_t + (t | expon(t,j)) * tensor(expon(t,j), {label_strs[i]}, {lie_group})
od
t = new_t
"""
        code += "print(t|z)\n"

    lines = lie_eval(code)
    try:
        return int(lines[0])
    except (IndexError, ValueError):
        return 0


def find_all_invariants(lie_group, dynkin_labels, max_degree=6):
    """Find which degree vectors have gauge-invariant structures (via tensor product).

    For each degree vector (n1, n2, ..., nk) with sum <= max_degree,
    checks if the raw tensor product R1^⊗n1 ⊗ ... ⊗ Rk^⊗nk contains a singlet.

    NOTE: This gives an upper bound — the raw tensor product doesn't account for
    identical-particle statistics. The actual PE count may be lower.
    Use this to identify WHICH degree vectors have primitives, then validate
    against PE for exact counts.

    Args:
        lie_group: LiE group name
        dynkin_labels: list of Dynkin label lists, one per rep type
        max_degree: max total degree to check

    Returns:
        dict mapping degree_tuple -> multiplicity of trivial in raw tensor product
    """
    from itertools import product as iprod
    k = len(dynkin_labels)
    label_strs = ['[' + ','.join(str(x) for x in dl) + ']' for dl in dynkin_labels]

    results = {}
    for degree in iprod(*(range(max_degree + 1) for _ in range(k))):
        total = sum(degree)
        if total < 2 or total > max_degree:
            continue

        rep_list = []
        for i, d in enumerate(degree):
            rep_list.extend([label_strs[i]] * d)

        if not rep_list:
            continue

        mult = singlet_in_multi_tensor(lie_group, *rep_list)
        if mult > 0:
            results[degree] = mult

    return results


if __name__ == '__main__':
    import sys

    # Test cases
    print("=== SU(2) fundamental [1] ===")
    invs = find_single_rep_invariants('A1', [1], max_degree=6)
    print(f"  Invariants: {invs}")
    print(f"  Contraction type: {get_contraction_type('A1', [1])}")

    print("\n=== SU(2) adjoint [2] ===")
    invs = find_single_rep_invariants('A1', [2], max_degree=6)
    print(f"  Invariants: {invs}")
    print(f"  Contraction type: {get_contraction_type('A1', [2])}")

    print("\n=== SU(3) fundamental [1,0] ===")
    invs = find_single_rep_invariants('A2', [1,0], max_degree=6)
    print(f"  Invariants: {invs}")
    # Should find: degree 3 alt (baryon = epsilon contraction)

    print("\n=== SU(3) adjoint [1,1] ===")
    invs = find_single_rep_invariants('A2', [1,1], max_degree=4)
    print(f"  Invariants: {invs}")

    print("\n=== G2 fundamental [1,0] ===")
    invs = find_single_rep_invariants('G2', [1,0], max_degree=6)
    print(f"  Invariants: {invs}")

    print("\n=== SO(5) = B2 vector [1,0] ===")
    invs = find_single_rep_invariants('B2', [1,0], max_degree=6)
    print(f"  Invariants: {invs}")
    print(f"  Contraction type: {get_contraction_type('B2', [1,0])}")

    print("\n=== Sp(4) = C2 fundamental [1,0] ===")
    invs = find_single_rep_invariants('C2', [1,0], max_degree=6)
    print(f"  Invariants: {invs}")
    print(f"  Contraction type: {get_contraction_type('C2', [1,0])}")

    print("\n=== Mixed: SU(3) fund [1,0] x antifund [0,1] ===")
    mixed = find_mixed_rep_invariants('A2', [[1,0], [0,1]])
    print(f"  Mixed invariants: {mixed}")

    print("\n=== Full catalog: SU(2) fund + adj ===")
    cat = build_invariant_catalog('A1', [[1], [2]], max_degree=6)
    print(f"  Single-rep: {cat['single_rep']}")
    print(f"  Mixed-rep: {cat['mixed_rep']}")
    print(f"  Contraction types: {cat['contraction_types']}")
