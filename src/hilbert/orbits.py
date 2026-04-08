"""
Orbit decomposition: group fields by equivalence under residual flavor symmetry.

Computes the automorphism group of the superpotential W to determine
which fields can be permuted while preserving W. Fields in the same
orbit are equivalent and yield the same deformed theory when used
as deformation operators.

Uses the bipartite graph structure from equivalence.py.
"""

import re
from collections import defaultdict
from itertools import permutations

import networkx as nx
from networkx.algorithms.isomorphism import GraphMatcher


def parse_fields_in_term(term_str):
    """Extract field names from a W term string like 'M01*r01*r03'."""
    return term_str.split('*')


def get_rep_type(field_name):
    """Extract rep type prefix from field name. 'r01' -> 'r0', 'M01' -> 'M0'."""
    if len(field_name) >= 2 and field_name[-1].isdigit():
        return field_name[:-1]
    return field_name


def build_w_graph(all_fields, w_terms):
    """Build a bipartite graph from the superpotential for automorphism analysis.

    Nodes: field nodes (labeled by rep type) + term nodes (labeled by structure)
    Edges: field appears in term (with multiplicity from powers)
    """
    G = nx.Graph()

    # Add field nodes with rep type label
    for f in all_fields:
        G.add_node(f"F_{f}", bipartite=0, rep_type=get_rep_type(f), field_name=f)

    # Add term nodes and edges
    for i, term in enumerate(w_terms):
        fields_in_term = parse_fields_in_term(term)
        term_node = f"T_{i}"
        # Term type = sorted tuple of rep types with multiplicities
        from collections import Counter
        rep_counts = tuple(sorted(Counter(get_rep_type(f) for f in fields_in_term).items()))
        G.add_node(term_node, bipartite=1, term_type=str(rep_counts))

        for f in fields_in_term:
            G.add_edge(term_node, f"F_{f}")

    return G


def find_field_automorphisms(all_fields, w_terms):
    """Find all field permutations that preserve W (automorphisms).

    Returns:
        list of dicts mapping field_name -> field_name (the permutation)
        Only permutations within the same rep type are considered.
    """
    if not w_terms:
        # W=0: all permutations within each rep type are automorphisms
        return None  # signal: full symmetry

    G = build_w_graph(all_fields, w_terms)

    def node_match(n1, n2):
        if n1.get('bipartite') != n2.get('bipartite'):
            return False
        if n1.get('bipartite') == 0:
            return n1.get('rep_type') == n2.get('rep_type')
        else:
            return n1.get('term_type') == n2.get('term_type')

    matcher = GraphMatcher(G, G, node_match=node_match)

    automorphisms = []
    for iso in matcher.isomorphisms_iter():
        # Extract the field-to-field mapping
        field_map = {}
        for k, v in iso.items():
            if k.startswith('F_') and v.startswith('F_'):
                f_from = k[2:]
                f_to = v[2:]
                if f_from != f_to:
                    field_map[f_from] = f_to
        if field_map:
            automorphisms.append(field_map)

    return automorphisms


def compute_orbits_from_automorphisms(all_fields, automorphisms, w_terms, singlet_prefixes=None):
    """Compute field orbits from the automorphism group.

    Two fields are in the same orbit if some automorphism maps one to the other.
    Uses union-find.
    """
    if singlet_prefixes is None:
        singlet_prefixes = set()

    # Union-Find
    parent = {f: f for f in all_fields}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(x, y):
        px, py = find(x), find(y)
        if px != py:
            parent[px] = py

    # Fields of the same rep type that aren't distinguished by W
    # are in the same orbit
    if automorphisms is None:
        # W=0: all copies of same rep type are one orbit
        by_type = defaultdict(list)
        for f in all_fields:
            by_type[get_rep_type(f)].append(f)
        for fields in by_type.values():
            for f in fields[1:]:
                union(fields[0], f)
    else:
        # Use automorphisms to determine equivalences
        # First: fields of same rep type that DON'T appear in any W term are equivalent
        w_fields = set()
        # (we need w_terms for this, but we only have automorphisms here)
        # Instead: two fields are equivalent if ANY automorphism maps one to the other
        for auto in automorphisms:
            for f_from, f_to in auto.items():
                if get_rep_type(f_from) == get_rep_type(f_to):
                    union(f_from, f_to)

        # Also: fields of same rep type not involved in any automorphism mapping
        # but that are also not distinguishable (spectators)
        # These would already be equivalent if there's an automorphism swapping them
        # If no automorphism involves them, they might still be equivalent
        # Actually, the identity automorphism leaves them fixed, but if field A
        # never maps to field B in any automorphism, they might still be in the
        # same orbit if there's a chain: A->C->B
        # The union-find handles transitivity.

        # Group true spectators: same rep type, not in any automorphism,
        # AND not appearing in any W term
        w_fields = set()
        for term in w_terms:
            for f in parse_fields_in_term(term):
                w_fields.add(f)

        by_type = defaultdict(list)
        for f in all_fields:
            by_type[get_rep_type(f)].append(f)
        for rt, fields in by_type.items():
            involved = set()
            for auto in automorphisms:
                for f_from, f_to in auto.items():
                    if get_rep_type(f_from) == rt:
                        involved.add(f_from)
                        involved.add(f_to)
            # True spectators: not involved in any automorphism AND not in W
            spectators = [f for f in fields if f not in involved and f not in w_fields]
            for f in spectators[1:]:
                union(spectators[0], f)

    # Build orbits from union-find
    orbit_map = defaultdict(list)
    for f in all_fields:
        orbit_map[find(f)].append(f)

    orbits = []
    for rep_field, fields in sorted(orbit_map.items()):
        rt = get_rep_type(fields[0])
        orbits.append({
            'fields': sorted(fields),
            'rep_type': rt,
            'is_singlet': rt in singlet_prefixes,
        })

    return orbits


def compute_orbits(all_fields, w_terms, singlet_prefixes=None):
    """Compute field orbits under residual flavor symmetry.

    Uses W automorphisms to determine which fields are equivalent.

    Args:
        all_fields: list of all field name strings
        w_terms: list of W term strings
        singlet_prefixes: set of field name prefixes that are singlets

    Returns:
        list of orbit dicts: {'fields': [...], 'rep_type': str, 'is_singlet': bool}
    """
    if singlet_prefixes is None:
        singlet_prefixes = set()

    if not w_terms:
        # W=0: all copies of same rep type form one orbit
        by_type = defaultdict(list)
        for f in all_fields:
            by_type[get_rep_type(f)].append(f)
        return [
            {'fields': sorted(fields), 'rep_type': rt, 'is_singlet': rt in singlet_prefixes}
            for rt, fields in sorted(by_type.items())
        ]

    # Compute automorphisms of W
    automorphisms = find_field_automorphisms(all_fields, w_terms)

    # Compute orbits from automorphisms
    return compute_orbits_from_automorphisms(all_fields, automorphisms, w_terms, singlet_prefixes)


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
    # Test cases
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

    # PE input test
    print("\n=== PE input for depth 1 ===")
    rep_map = {'r0': [1]}
    labels, mults, oidx = orbits_to_pe_input(orbits1, rep_map)
    print(f"  Dynkin labels: {labels}")
    print(f"  Multiplicities: {mults}")
    print(f"  PE partitions: (7)^{len(labels)} = {7**len(labels)}")
