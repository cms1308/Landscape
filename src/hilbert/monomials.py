"""
Expand PE degree vectors into explicit gauge-invariant monomials.

Given:
- A PE degree vector (n_1, n_2, ..., n_k) where n_i = number of fields from orbit i
- Orbit structure: for each orbit, the list of field names and the contraction type
- Invariant catalog from invariants.py

Produces:
- All explicit monomial strings (e.g. "r01*r03", "M01*r03*r04")
"""

from itertools import combinations, combinations_with_replacement, product as iprod


def select_from_orbit(orbit_fields, count, contraction_type):
    """Generate all ways to select `count` fields from an orbit.

    Args:
        orbit_fields: list of field name strings (e.g. ["r03","r04","r05","r06","r07","r08"])
        count: number of fields to select
        contraction_type: "alt" (antisymmetric), "sym" (symmetric), or "any" (all combos)

    Returns:
        list of tuples of field names
    """
    if count == 0:
        return [()]
    if count > len(orbit_fields) and contraction_type == "alt":
        return []  # can't antisymmetrize more fields than available

    if contraction_type == "alt":
        return list(combinations(orbit_fields, count))
    elif contraction_type == "sym":
        return list(combinations_with_replacement(orbit_fields, count))
    else:  # "any" — all ordered selections (for mixed-rep products)
        return list(combinations_with_replacement(orbit_fields, count))


def build_monomial(field_tuple):
    """Convert a tuple of field names to a monomial string."""
    if not field_tuple:
        return ""
    return "*".join(field_tuple)


def expand_single_degree(orbits, degree_vector, invariant_catalog):
    """Expand a single PE degree vector into explicit monomials.

    Args:
        orbits: list of dicts with keys:
            'fields': list of field name strings
            'rep_idx': index into the invariant catalog's rep list
            'contraction': "alt", "sym", or "any"
            'is_singlet': bool
        degree_vector: tuple of ints (n_1, n_2, ...) — how many fields from each orbit
        invariant_catalog: from invariants.build_invariant_catalog

    Returns:
        list of monomial strings
    """
    # For each orbit, generate selections
    per_orbit_selections = []
    for i, orbit in enumerate(orbits):
        d = degree_vector[i]
        if d == 0:
            per_orbit_selections.append([()])
            continue

        if orbit.get('is_singlet', False):
            # Singlet: only one field, power d
            if d == 1:
                per_orbit_selections.append([(orbit['fields'][0],)])
            else:
                # d copies of the same singlet field: field^d
                per_orbit_selections.append([tuple([orbit['fields'][0]] * d)])
        else:
            ct = orbit['contraction']
            selections = select_from_orbit(orbit['fields'], d, ct)
            if not selections:
                return []  # impossible degree for this orbit
            per_orbit_selections.append(selections)

    # Cartesian product across orbits
    monomials = []
    for combo in iprod(*per_orbit_selections):
        fields = []
        for selection in combo:
            fields.extend(selection)
        if fields:
            monomials.append(build_monomial(tuple(fields)))

    return monomials


def expand_all_degrees(orbits, pe_degrees, invariant_catalog, r_charges, max_R=2.0):
    """Expand all PE degree vectors into explicit monomials, filtered by R-charge.

    Args:
        orbits: list of orbit dicts (see expand_single_degree)
        pe_degrees: list of (degree_tuple, multiplicity) from PE
        invariant_catalog: from invariants.build_invariant_catalog
        r_charges: dict mapping field name -> R-charge (float)
        max_R: only include operators with total R < max_R

    Returns:
        list of dicts with 'monomial', 'R', 'degree'
    """
    all_operators = []
    seen = set()

    for degree, pe_mult in pe_degrees:
        if all(d == 0 for d in degree):
            continue

        # Quick R-charge check: compute min possible R for this degree
        # (using min R per orbit)
        min_R = 0
        for i, orbit in enumerate(orbits):
            if degree[i] > 0 and orbit['fields']:
                min_field_R = min(r_charges.get(f, 0) for f in orbit['fields'])
                min_R += degree[i] * min_field_R
        if min_R >= max_R:
            continue

        monomials = expand_single_degree(orbits, degree, invariant_catalog)

        for mono in monomials:
            if mono in seen:
                continue
            seen.add(mono)

            # Compute R-charge
            fields = mono.split('*')
            R = sum(r_charges.get(f, 0) for f in fields)
            if R < max_R:
                all_operators.append({
                    'monomial': mono,
                    'R': R,
                    'degree': degree,
                })

    return sorted(all_operators, key=lambda x: x['R'])


def build_product_monomials(base_monomials, r_charges, max_R=2.0, max_products=2):
    """Build products of base monomials (e.g. meson*meson = quartic).

    Args:
        base_monomials: list of monomial strings (the "generators")
        r_charges: dict field -> R-charge
        max_R: filter
        max_products: max number of factors in a product (2 = quartic from mesons)

    Returns:
        list of dicts with 'monomial', 'R'
    """
    products = []
    seen = set()

    for n_factors in range(2, max_products + 1):
        for combo in combinations_with_replacement(range(len(base_monomials)), n_factors):
            fields = []
            for idx in combo:
                fields.extend(base_monomials[idx].split('*'))
            fields_sorted = tuple(sorted(fields))
            if fields_sorted in seen:
                continue

            R = sum(r_charges.get(f, 0) for f in fields)
            if R >= max_R:
                continue

            seen.add(fields_sorted)
            mono = '*'.join(fields)
            products.append({'monomial': mono, 'R': R})

    return sorted(products, key=lambda x: x['R'])


def generate_all_operators(orbits, pe_degrees, invariant_catalog, r_charges,
                           singlet_fields=None, max_R=2.0):
    """Full operator generation: PE monomials + singlet products + quartic products.

    Args:
        orbits: orbit structure for gauge-charged fields
        pe_degrees: from Hilbert series PE
        invariant_catalog: from invariants.py
        r_charges: dict field -> R-charge (all fields including singlets)
        singlet_fields: list of singlet field names (e.g. ["M01"])
        max_R: R-charge cutoff

    Returns:
        list of operator dicts {'monomial', 'R'}, sorted by R
    """
    if singlet_fields is None:
        singlet_fields = []

    # 1. Expand PE degrees into gauge-invariant monomials
    gauge_ops = expand_all_degrees(orbits, pe_degrees, invariant_catalog, r_charges, max_R)
    gauge_mono_strs = [op['monomial'] for op in gauge_ops]

    # 2. Singlet fields themselves
    singlet_ops = []
    for sf in singlet_fields:
        R = r_charges.get(sf, 0)
        if R < max_R:
            singlet_ops.append({'monomial': sf, 'R': R})
        # Singlet powers
        for power in range(2, int(max_R / R) + 1 if R > 0 else 2):
            Rp = power * R
            if Rp < max_R:
                singlet_ops.append({'monomial': '*'.join([sf] * power), 'R': Rp})

    # 3. Singlet * gauge-invariant products
    mixed_ops = []
    for sf in singlet_fields:
        R_sf = r_charges.get(sf, 0)
        for gop in gauge_ops:
            R_total = R_sf + gop['R']
            if R_total < max_R:
                mixed_ops.append({
                    'monomial': sf + '*' + gop['monomial'],
                    'R': R_total,
                })
        # singlet^2 * gauge ops
        for gop in gauge_ops:
            R_total = 2 * R_sf + gop['R']
            if R_total < max_R:
                mixed_ops.append({
                    'monomial': sf + '*' + sf + '*' + gop['monomial'],
                    'R': R_total,
                })

    # 4. Products of gauge-invariant monomials (quartic, etc.)
    product_ops = build_product_monomials(gauge_mono_strs, r_charges, max_R, max_products=2)

    # 5. Singlet * product ops
    singlet_product_ops = []
    for sf in singlet_fields:
        R_sf = r_charges.get(sf, 0)
        for pop in product_ops:
            R_total = R_sf + pop['R']
            if R_total < max_R:
                singlet_product_ops.append({
                    'monomial': sf + '*' + pop['monomial'],
                    'R': R_total,
                })

    # Combine all, deduplicate
    all_ops = gauge_ops + singlet_ops + mixed_ops + product_ops + singlet_product_ops
    seen = set()
    unique_ops = []
    for op in all_ops:
        key = tuple(sorted(op['monomial'].split('*')))
        if key not in seen:
            seen.add(key)
            unique_ops.append(op)

    return sorted(unique_ops, key=lambda x: x['R'])


if __name__ == '__main__':
    # Test: SU(2) 8 fund at depth 1, W = M01*r01*r02
    # Orbits: {r01}, {r02}, {r03..r08}, {M01}
    # R-charges: r01=r02=0.5403, r03-08=0.4866, M01=0.9194

    r_charges = {
        'r01': 0.5403, 'r02': 0.5403,
        'r03': 0.4866, 'r04': 0.4866, 'r05': 0.4866,
        'r06': 0.4866, 'r07': 0.4866, 'r08': 0.4866,
        'M01': 0.9194,
    }

    # Orbits for gauge-charged fields (fund of SU(2) = pseudo-real, contraction = alt)
    orbits = [
        {'fields': ['r01'], 'rep_idx': 0, 'contraction': 'alt', 'is_singlet': False},
        {'fields': ['r02'], 'rep_idx': 0, 'contraction': 'alt', 'is_singlet': False},
        {'fields': ['r03', 'r04', 'r05', 'r06', 'r07', 'r08'],
         'rep_idx': 0, 'contraction': 'alt', 'is_singlet': False},
    ]

    # PE degrees: from per-type PE with 3 orbit types
    # For SU(2) fund^2: meson from any two orbits
    # Simulate: the PE would give degrees like (1,1,0), (1,0,1), (0,1,1), (0,0,2), etc.
    # with multiplicities
    pe_degrees = [
        ((1, 1, 0), 1),   # r01*r02 — but this is killed by F-term!
        ((1, 0, 1), 1),   # r01*r0i (6 choices from orbit 2)
        ((0, 1, 1), 1),   # r02*r0i (6 choices)
        ((0, 0, 2), 1),   # r0i*r0j from spectator orbit (C(6,2)=15 choices)
        ((2, 0, 2), 1),   # quartic involving r01
        ((0, 2, 2), 1),   # quartic involving r02
        ((1, 1, 2), 1),   # quartic involving both r01, r02
        ((0, 0, 4), 1),   # pure spectator quartic
    ]

    ops = generate_all_operators(
        orbits, pe_degrees, {}, r_charges,
        singlet_fields=['M01'], max_R=2.0
    )

    print(f"Total operators with R < 2: {len(ops)}")
    print("\nBy R-charge range:")
    ranges = {}
    for op in ops:
        bucket = f"{op['R']:.2f}"
        ranges.setdefault(bucket, []).append(op['monomial'])

    for r_val in sorted(ranges.keys()):
        print(f"  R~{r_val}: {len(ranges[r_val])} ops")
        for mono in ranges[r_val][:3]:
            print(f"    {mono}")
        if len(ranges[r_val]) > 3:
            print(f"    ... ({len(ranges[r_val])-3} more)")
