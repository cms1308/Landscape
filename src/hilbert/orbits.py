"""
Orbit decomposition via graph isomorphism (matching Nf=2N.nb approach).

For each field, compute the superpotential graph with that field appended
as a marker. Fields that produce isomorphic graphs are in the same orbit
(equivalent under residual flavor symmetry).

This is the toGlobalGraph approach from Nf=2N.nb:
  pairs = GroupBy[funds, toGlobalGraph[Append[w, ToString[#]]] &]
"""

from collections import defaultdict
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from deform.equivalence import canonical_hash, build_theory_graph, are_equivalent


def get_rep_type(field_name):
    """Extract rep type prefix. 'r01' -> 'r0', 'M01' -> 'M0'."""
    if len(field_name) >= 2 and field_name[-1].isdigit():
        return field_name[:-1]
    return field_name


def compute_orbits(all_fields, w_terms, singlet_prefixes=None):
    """Compute field orbits via graph isomorphism.

    For each field f (of the same rep type), compute the graph of W ∪ {f}
    and group fields with isomorphic graphs.

    Args:
        all_fields: list of all field name strings
        w_terms: list of W term strings
        singlet_prefixes: set of field name prefixes that are singlets

    Returns:
        list of orbit dicts: {'fields': [...], 'rep_type': str, 'is_singlet': bool}
    """
    if singlet_prefixes is None:
        singlet_prefixes = set()

    # Group fields by rep type first
    by_type = defaultdict(list)
    for f in all_fields:
        by_type[get_rep_type(f)].append(f)

    orbits = []

    for rep_type, fields in sorted(by_type.items()):
        is_singlet = rep_type in singlet_prefixes

        if not w_terms:
            # W=0: all copies equivalent
            orbits.append({
                'fields': sorted(fields),
                'rep_type': rep_type,
                'is_singlet': is_singlet,
            })
            continue

        # For each field, compute canonical hash of W ∪ {field_name}
        # Fields with the same hash are candidates for the same orbit
        hash_groups = defaultdict(list)
        for f in fields:
            w_extended = w_terms + [f]
            h = canonical_hash(w_extended)
            hash_groups[h].append(f)

        # Within each hash group, verify by graph isomorphism (hash may have collisions)
        for h, group in hash_groups.items():
            sub_orbits = []
            for f in group:
                placed = False
                for so in sub_orbits:
                    rep_f = so[0]
                    w_ext_f = w_terms + [f]
                    w_ext_rep = w_terms + [rep_f]
                    if are_equivalent(w_ext_f, w_ext_rep):
                        so.append(f)
                        placed = True
                        break
                if not placed:
                    sub_orbits.append([f])

            for so in sub_orbits:
                orbits.append({
                    'fields': sorted(so),
                    'rep_type': rep_type,
                    'is_singlet': is_singlet,
                })

    return orbits


def orbits_to_pe_input(orbits, rep_type_to_dynkin):
    """Convert orbits to input for compute_hilbert_series.

    Returns:
        (dynkin_labels, multiplicities, orbit_indices)
    """
    dynkin_labels = []
    multiplicities = []
    orbit_indices = []

    for i, orbit in enumerate(orbits):
        if orbit['is_singlet']:
            continue
        dl = rep_type_to_dynkin.get(orbit['rep_type'])
        if dl is None:
            continue
        dynkin_labels.append(dl)
        multiplicities.append(len(orbit['fields']))
        orbit_indices.append(i)

    return dynkin_labels, multiplicities, orbit_indices


if __name__ == '__main__':
    import sys, os
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

    all_fields = [f'r0{i}' for i in range(1, 9)] + ['M01']

    print("=== W = {} (depth 0) ===")
    orbits0 = compute_orbits(all_fields, [], singlet_prefixes={'M0'})
    for o in orbits0:
        print(f"  {o['rep_type']}: {o['fields']} (singlet={o['is_singlet']})")

    print("\n=== W = {M01*r01*r02} (depth 1) ===")
    orbits1 = compute_orbits(all_fields, ['M01*r01*r02'], singlet_prefixes={'M0'})
    for o in orbits1:
        print(f"  {o['rep_type']}: {o['fields']} (singlet={o['is_singlet']})")

    print("\n=== W = {M01*r01*r02, r03*r04} (depth 2) ===")
    orbits2 = compute_orbits(all_fields, ['M01*r01*r02', 'r03*r04'], singlet_prefixes={'M0'})
    for o in orbits2:
        print(f"  {o['rep_type']}: {o['fields']} (singlet={o['is_singlet']})")

    print("\n=== W = {M01*r01*r02, r01*r03} (depth 2, asymmetric) ===")
    orbits3 = compute_orbits(all_fields, ['M01*r01*r02', 'r01*r03'], singlet_prefixes={'M0'})
    for o in orbits3:
        print(f"  {o['rep_type']}: {o['fields']} (singlet={o['is_singlet']})")

    print("\n=== PE input for depth 1 ===")
    rep_map = {'r0': [1]}
    labels, mults, oidx = orbits_to_pe_input(orbits1, rep_map)
    print(f"  Dynkin labels: {labels}")
    print(f"  Multiplicities: {mults}")
    print(f"  PE partitions: (7)^{len(labels)} = {7**len(labels)}")
