"""
Theory equivalence checking via superpotential graph isomorphism.

Port of toGlobalGraph from Nf=2N.nb.

Two theories (same gauge group) are equivalent if there exists a relabeling
of fields (within each representation type) that maps one superpotential
to the other. This is checked via bipartite graph isomorphism.

Graph construction:
- Term nodes: one per superpotential term, labeled by type "T"
- Field nodes: one per field variable, labeled by rep type (e.g. "q", "phi", "M")
- Edges: connect term node to field node if the field appears in that term
  Edge multiplicity = power of the field in the term
- Two theories are equivalent iff their graphs are isomorphic with
  node labels (types) preserved.
"""

import re
from collections import Counter
import networkx as nx
from networkx.algorithms.isomorphism import GraphMatcher, MultiGraphMatcher


def parse_monomial(term_str):
    """Parse a monomial string like 'r01*r02*r11' or 'r06^2*r07^2' into {field: power} dict."""
    result = Counter()
    for factor in term_str.split('*'):
        if '^' in factor:
            base, exp = factor.split('^', 1)
            result[base] += int(exp)
        else:
            result[factor] += 1
    return dict(result)


def extract_rep_type(field_name):
    """Extract the representation type from a field name.

    Convention: field = repname + single_digit_copy_index
    e.g. 'r01' -> 'r0', 'r11' -> 'r1', 'M21' -> 'M2', 'phi1' -> 'phi'
    The copy index is the LAST digit only.
    """
    if len(field_name) >= 2 and field_name[-1].isdigit():
        return field_name[:-1]
    return field_name


def build_theory_graph(w_terms, matter_info=None):
    """Build a bipartite graph from superpotential terms.

    Args:
        w_terms: list of monomial strings, e.g. ["r01*r02", "M1*r01*r11"]
        matter_info: optional dict mapping field names to rep type labels

    Returns:
        networkx.Graph with node attributes for matching
    """
    G = nx.MultiGraph()

    # Collect all fields
    all_fields = set()
    parsed_terms = []
    for term in w_terms:
        parsed = parse_monomial(term)
        parsed_terms.append(parsed)
        all_fields.update(parsed.keys())

    # Add field nodes with rep type attribute
    for field in sorted(all_fields):
        rep_type = extract_rep_type(field)
        G.add_node(f"F_{field}", bipartite=0, node_type=rep_type)

    # Add term nodes
    for i, parsed in enumerate(parsed_terms):
        term_label = f"T_{i}"
        # Term type = sorted tuple of (rep_type, power) — this captures the structure
        term_type = tuple(sorted(
            (extract_rep_type(f), p) for f, p in parsed.items()
        ))
        G.add_node(term_label, bipartite=1, node_type=str(term_type))

        # Add edges with multiplicity
        for field, power in parsed.items():
            for _ in range(power):
                G.add_edge(term_label, f"F_{field}")

    return G


def node_match(n1_attrs, n2_attrs):
    """Node matching function: types must agree."""
    return n1_attrs.get('node_type') == n2_attrs.get('node_type')


def are_equivalent(w1, w2, matter1=None, matter2=None):
    """Check if two theories (same gauge group, same matter multiplicities) are equivalent.

    Args:
        w1, w2: lists of superpotential term strings
        matter1, matter2: optional matter info dicts

    Returns:
        bool: True if theories are equivalent up to field relabeling
    """
    if len(w1) != len(w2):
        return False

    if not w1 and not w2:
        return True  # both W=0

    G1 = build_theory_graph(w1, matter1)
    G2 = build_theory_graph(w2, matter2)

    # Quick check: same number of nodes and edges
    if G1.number_of_nodes() != G2.number_of_nodes():
        return False
    if G1.number_of_edges() != G2.number_of_edges():
        return False

    # Graph isomorphism with node type matching
    matcher = MultiGraphMatcher(G1, G2, node_match=node_match)
    return matcher.is_isomorphic()


def canonical_hash(w_terms):
    """Compute a hash for quick equivalence pre-check.

    Two equivalent theories must have the same hash (but not vice versa).
    """
    if not w_terms:
        return (0,)

    term_signatures = []
    for term in w_terms:
        parsed = parse_monomial(term)
        # Signature: sorted tuple of (rep_type, power)
        sig = tuple(sorted(
            (extract_rep_type(f), p) for f, p in parsed.items()
        ))
        term_signatures.append(sig)

    return tuple(sorted(term_signatures))


def deduplicate_theories(theories):
    """Remove duplicate theories from a list.

    Args:
        theories: list of theory dicts, each with 'w' key (list of W terms)

    Returns:
        list of unique theories
    """
    if not theories:
        return []

    # Group by canonical hash for fast pre-filtering
    by_hash = {}
    for t in theories:
        h = canonical_hash(t.get('w', []))
        if h not in by_hash:
            by_hash[h] = []
        by_hash[h].append(t)

    unique = []
    for h, group in by_hash.items():
        # Within each hash group, check pairwise equivalence
        representatives = []
        for t in group:
            is_dup = False
            for rep in representatives:
                if are_equivalent(t.get('w', []), rep.get('w', [])):
                    is_dup = True
                    break
            if not is_dup:
                representatives.append(t)
        unique.extend(representatives)

    return unique


def deduplicate_operators(operators, w_terms):
    """Deduplicate operators by comparing toGlobalGraph(W + [op]).

    Two operators are equivalent if appending them to W gives isomorphic graphs.
    This is the Nf=2N.nb approach:
      DeleteDuplicatesBy[ops, toGlobalGraph[Append[w, ToString[#]]] &]

    Args:
        operators: list of operator dicts with 'monomial' key
        w_terms: current W terms (list of strings)

    Returns:
        list of unique operator dicts (one representative per equivalence class)
    """
    if not operators:
        return []

    # Group by canonical hash of W + [op]
    by_hash = {}
    for op in operators:
        w_ext = w_terms + [op['monomial']]
        h = canonical_hash(w_ext)
        if h not in by_hash:
            by_hash[h] = []
        by_hash[h].append(op)

    unique = []
    for h, group in by_hash.items():
        representatives = []
        for op in group:
            is_dup = False
            w_ext = w_terms + [op['monomial']]
            for rep in representatives:
                w_ext_rep = w_terms + [rep['monomial']]
                if are_equivalent(w_ext, w_ext_rep):
                    is_dup = True
                    break
            if not is_dup:
                representatives.append(op)
        unique.extend(representatives)

    return unique


if __name__ == '__main__':
    # Test cases
    print("=== Test: equivalence checking ===")

    # Same theory, different field labeling
    w1 = ["r01*r02", "M1*r01*r03"]
    w2 = ["r01*r03", "M1*r01*r02"]  # relabel r02 <-> r03
    print(f"w1={w1}")
    print(f"w2={w2}")
    print(f"Equivalent? {are_equivalent(w1, w2)}")  # Should be True

    # Different theory
    w3 = ["r01*r02", "M1*r02*r03"]
    print(f"\nw3={w3}")
    print(f"w1 ~ w3? {are_equivalent(w1, w3)}")  # Should be False (different structure)

    # Same hash but different
    w4 = ["r01*r02", "r03*r04"]
    w5 = ["r01*r03", "r02*r04"]
    print(f"\nw4={w4}")
    print(f"w5={w5}")
    print(f"hash(w4)={canonical_hash(w4)}")
    print(f"hash(w5)={canonical_hash(w5)}")
    print(f"w4 ~ w5? {are_equivalent(w4, w5)}")  # Should be True (relabel)

    # Deduplication
    print("\n=== Test: deduplication ===")
    theories = [
        {'w': ["r01*r02"], 'desc': 'A'},
        {'w': ["r01*r03"], 'desc': 'B'},  # same as A up to relabeling
        {'w': ["r01*r11"], 'desc': 'C'},  # different (different rep types)
    ]
    unique = deduplicate_theories(theories)
    print(f"Input: {len(theories)} theories")
    print(f"Unique: {len(unique)} theories")
    for t in unique:
        print(f"  {t['desc']}: {t['w']}")
