"""
Orbit decomposition: group fields by equivalence under residual flavor symmetry.

For W=0: all copies of the same rep type form one orbit.
For W≠0: fields appearing in W are distinguished (singleton orbits);
remaining copies of the same type form one orbit.
"""

import re
from collections import defaultdict


def parse_fields_in_term(term_str):
    """Extract field names from a W term string like 'M01*r01*r03'."""
    return term_str.split('*')


def get_rep_type(field_name):
    """Extract rep type prefix from field name. E.g. 'r01' -> 'r0', 'M01' -> 'M0'."""
    if len(field_name) >= 2 and field_name[-1].isdigit():
        return field_name[:-1]
    return field_name


def compute_orbits(all_fields, w_terms, singlet_prefixes=None):
    """Compute field orbits under residual flavor symmetry.

    Args:
        all_fields: list of all field name strings (e.g. ['r01','r02',...,'r08','M01'])
        w_terms: list of W term strings (e.g. ['M01*r01*r02'])
        singlet_prefixes: set of field name prefixes that are singlets (e.g. {'M0','X0'})

    Returns:
        list of orbit dicts:
            'fields': list of field names in this orbit
            'rep_type': rep type prefix
            'is_singlet': bool
    """
    if singlet_prefixes is None:
        singlet_prefixes = set()

    # Find all fields appearing in W terms
    w_fields = set()
    for term in w_terms:
        for f in parse_fields_in_term(term):
            w_fields.add(f)

    # Group fields by rep type
    by_type = defaultdict(list)
    for f in all_fields:
        by_type[get_rep_type(f)].append(f)

    orbits = []

    for rep_type, fields in sorted(by_type.items()):
        is_singlet = rep_type in singlet_prefixes

        if not w_terms:
            # W=0: all copies are one orbit
            orbits.append({
                'fields': fields,
                'rep_type': rep_type,
                'is_singlet': is_singlet,
            })
        else:
            # W≠0: split into distinguished (in W) and spectators
            distinguished = [f for f in fields if f in w_fields]
            spectators = [f for f in fields if f not in w_fields]

            # Each distinguished field is its own orbit
            for f in distinguished:
                orbits.append({
                    'fields': [f],
                    'rep_type': rep_type,
                    'is_singlet': is_singlet,
                })

            # All spectators of same type form one orbit
            if spectators:
                orbits.append({
                    'fields': spectators,
                    'rep_type': rep_type,
                    'is_singlet': is_singlet,
                })

    return orbits


def orbits_to_pe_input(orbits, rep_type_to_dynkin):
    """Convert orbits to input for compute_hilbert_series.

    Args:
        orbits: list of orbit dicts (from compute_orbits)
        rep_type_to_dynkin: dict mapping rep_type prefix to Dynkin label list

    Returns:
        (dynkin_labels, multiplicities) for non-singlet orbits only
        Also returns orbit_indices: mapping from PE type index to orbit index
    """
    dynkin_labels = []
    multiplicities = []
    orbit_indices = []  # which orbit each PE type corresponds to

    for i, orbit in enumerate(orbits):
        if orbit['is_singlet']:
            continue  # singlets don't enter gauge projection
        dl = rep_type_to_dynkin.get(orbit['rep_type'])
        if dl is None:
            continue
        dynkin_labels.append(dl)
        multiplicities.append(len(orbit['fields']))
        orbit_indices.append(i)

    return dynkin_labels, multiplicities, orbit_indices


if __name__ == '__main__':
    # Test: SU(2) 8 fund + 1 singlet M, W = M01*r01*r02
    all_fields = [f'r0{i}' for i in range(1, 9)] + ['M01']

    print("=== W = {} (depth 0) ===")
    orbits0 = compute_orbits(all_fields, [], singlet_prefixes={'M0'})
    for o in orbits0:
        print(f"  {o['rep_type']}: {o['fields']} (singlet={o['is_singlet']})")

    print("\n=== W = {M01*r01*r02} (depth 1) ===")
    orbits1 = compute_orbits(all_fields, ['M01*r01*r02'], singlet_prefixes={'M0'})
    for o in orbits1:
        print(f"  {o['rep_type']}: {o['fields']} (singlet={o['is_singlet']})")

    print("\n=== PE input for depth 1 ===")
    rep_map = {'r0': [1]}  # SU(2) fundamental
    labels, mults, oidx = orbits_to_pe_input(orbits1, rep_map)
    print(f"  Dynkin labels: {labels}")
    print(f"  Multiplicities: {mults}")
    print(f"  Orbit indices: {oidx}")
    print(f"  PE partitions: (order+1)^{len(labels)} = 7^{len(labels)} = {7**len(labels)}")
