"""
GIO Cache: stores primitives, products, and PE validation for each (group, matter).

Per (gauge_group, rep_types, multiplicities):
  - primitives: list of (degree, sym_type, monomials)
  - products: dict of order -> list of monomials at that order
  - pe_counts: dict of order -> {degree: PE_count}
  - matching: dict of order -> {degree: (our_count, pe_count, status)}

The cache is built once per seed theory and extended when higher orders are needed.
"""

import os
import json
import math
import ast
from collections import defaultdict, Counter
from itertools import combinations, combinations_with_replacement

import sys
_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _dir)
sys.path.insert(0, os.path.join(_dir, '..'))

from invariants import find_single_rep_invariants, find_mixed_rep_invariants
from hilbert import compute_hilbert_series


class GIOCache:
    """Cache of gauge-invariant operators for a specific (group, matter) combination."""

    def __init__(self, lie_group, rep_fields, dynkin_labels, singlet_fields=None):
        """
        Args:
            lie_group: LiE group name (e.g. "A1", "B3")
            rep_fields: dict {rep_name: [field_name, ...]} for gauge-charged reps
            dynkin_labels: list of Dynkin label lists, one per rep name (same order as rep_fields keys)
            singlet_fields: list of singlet field names
        """
        self.lie_group = lie_group
        self.rep_names = list(rep_fields.keys())
        self.rep_fields = rep_fields
        self.dynkin_labels = dynkin_labels
        self.singlet_fields = singlet_fields or []

        # All gauge fields flat
        self.all_gauge_fields = []
        for name in self.rep_names:
            self.all_gauge_fields.extend(rep_fields[name])

        # Per-rep PE input
        self.pe_labels = []
        self.pe_mults = []
        self._label_to_idx = {}
        for i, name in enumerate(self.rep_names):
            dl = tuple(dynkin_labels[i])
            if dl not in self._label_to_idx:
                self._label_to_idx[dl] = len(self.pe_labels)
                self.pe_labels.append(list(dl))
                self.pe_mults.append(0)
            self.pe_mults[self._label_to_idx[dl]] += len(rep_fields[name])

        # Cache storage
        self.primitives = None  # list of (degree, sym_type, [monomial_strs])
        self.products_by_order = {}  # order -> list of monomial strs (gauge-only)
        self.pe_counts_by_order = {}  # order -> {degree_tuple: count}
        self.matching_by_order = {}  # order -> {degree_tuple: (ours, pe, status)}
        self.max_order_built = 0

    def _find_primitives(self):
        """Find primitive gauge invariants via invariants.py."""
        self.primitives = []

        # Single-rep primitives
        for i, name in enumerate(self.rep_names):
            dl = self.dynkin_labels[i]
            invs = find_single_rep_invariants(self.lie_group, dl, max_degree=8)
            for degree, sym_type, mult in invs:
                if sym_type == "alt":
                    monos = ['*'.join(c) for c in combinations(self.rep_fields[name], degree)]
                else:
                    monos = ['*'.join(c) for c in combinations_with_replacement(self.rep_fields[name], degree)]
                if monos:
                    self.primitives.append((degree, sym_type, name, monos))

        # Mixed-rep primitives (bilinear between different rep types)
        mixed = find_mixed_rep_invariants(self.lie_group, self.dynkin_labels)
        for i, j, mult in mixed:
            if i == j:
                continue
            monos = [f"{fi}*{fj}"
                     for fi in self.rep_fields[self.rep_names[i]]
                     for fj in self.rep_fields[self.rep_names[j]]]
            if monos:
                self.primitives.append((2, "mixed", f"{self.rep_names[i]}x{self.rep_names[j]}", monos))

    def _get_all_primitive_monomials(self):
        """Return flat list of all primitive monomial strings."""
        if self.primitives is None:
            self._find_primitives()
        monos = []
        for _, _, _, mono_list in self.primitives:
            monos.extend(mono_list)
        return monos

    def _build_products(self, max_degree):
        """Build products of primitives up to max_degree (in total field count)."""
        prims = self._get_all_primitive_monomials()

        # Track by sorted field tuple to avoid duplicates
        seen = set()
        products = []

        # Primitives themselves
        for p in prims:
            fields = p.split('*')
            if len(fields) <= max_degree:
                key = tuple(sorted(fields))
                if key not in seen:
                    seen.add(key)
                    products.append(p)

        # Products of 2 primitives
        for i in range(len(prims)):
            fi = prims[i].split('*')
            if len(fi) * 2 > max_degree:
                continue
            for j in range(i, len(prims)):
                fj = prims[j].split('*')
                total_deg = len(fi) + len(fj)
                if total_deg > max_degree:
                    continue
                combined = fi + fj
                key = tuple(sorted(combined))
                if key not in seen:
                    seen.add(key)
                    products.append('*'.join(combined))

        # Products of 3 primitives
        if max_degree >= 6:
            for i in range(len(prims)):
                fi = prims[i].split('*')
                if len(fi) * 3 > max_degree:
                    continue
                for j in range(i, len(prims)):
                    fj = prims[j].split('*')
                    if len(fi) + len(fj) * 2 > max_degree:
                        continue
                    for k in range(j, len(prims)):
                        fk = prims[k].split('*')
                        total_deg = len(fi) + len(fj) + len(fk)
                        if total_deg > max_degree:
                            continue
                        combined = fi + fj + fk
                        key = tuple(sorted(combined))
                        if key not in seen:
                            seen.add(key)
                            products.append('*'.join(combined))

        return products

    def _compute_pe(self, order):
        """Compute per-rep PE at given order."""
        if not self.pe_labels:
            return {}
        hs = compute_hilbert_series(self.lie_group, self.pe_labels, self.pe_mults,
                                     order=order, compute_pl=False)
        pe_degrees = {}
        for item_str in hs['pe_raw']:
            matrix = ast.literal_eval(item_str)
            rows = matrix if isinstance(matrix[0], list) else [matrix]
            for row in rows:
                degree = tuple(row[:-1])
                mult = row[-1]
                if any(d > 0 for d in degree):
                    pe_degrees[degree] = mult
        return pe_degrees

    def _count_products_by_degree(self, products):
        """Count gauge-only products by per-rep degree vector."""
        counts = defaultdict(int)
        for op in products:
            fields = op.split('*')
            degree = [0] * len(self.pe_labels)
            for f in fields:
                for name in self.rep_names:
                    if f in self.rep_fields[name]:
                        dl = tuple(self.dynkin_labels[self.rep_names.index(name)])
                        degree[self._label_to_idx[dl]] += 1
                        break
            counts[tuple(degree)] += 1
        return dict(counts)

    def build_and_validate(self, order):
        """Build products up to `order` and validate against PE.

        Returns:
            True if validation passed, False otherwise.
        """
        if self.primitives is None:
            self._find_primitives()

        # Build products
        products = self._build_products(max_degree=order)
        self.products_by_order[order] = products

        # Compute PE
        pe_counts = self._compute_pe(order)
        self.pe_counts_by_order[order] = pe_counts

        # Count our products by degree
        our_counts = self._count_products_by_degree(products)

        # Validate: exact match (our <= PE, our > 0 where PE > 0)
        matching = {}
        passed = True
        for degree, pe_mult in pe_counts.items():
            our_mult = our_counts.get(degree, 0)
            if our_mult > pe_mult:
                status = "EXTRA"
                passed = False
            elif our_mult == 0 and pe_mult > 0:
                status = "MISSING"
                passed = False
            elif our_mult == pe_mult:
                status = "EXACT"
            else:
                status = "OK"  # our < PE from multiplicity
            matching[str(degree)] = (our_mult, pe_mult, status)

        # Check for extras not in PE
        for degree, our_mult in our_counts.items():
            if degree not in pe_counts and our_mult > 0:
                matching[str(degree)] = (our_mult, 0, "EXTRA")
                passed = False

        self.matching_by_order[order] = matching
        self.max_order_built = max(self.max_order_built, order)

        if passed:
            print(f"  GIO cache: order {order}, {len(products)} products, PE validation PASSED")
        else:
            print(f"  GIO cache: order {order}, {len(products)} products, PE validation FAILED")
            for deg_str, (ours, pe, status) in matching.items():
                if status in ("MISSING", "EXTRA"):
                    print(f"    {status} at degree {deg_str}: ours={ours}, PE={pe}")

        return passed

    def get_operators(self, r_charges, max_R=2.0, w_terms=None):
        """Get all GIO monomials (gauge + singlet products) filtered by R < max_R.

        Automatically extends order if needed based on min R-charge.

        Args:
            r_charges: dict field -> R-charge
            max_R: R-charge cutoff
            w_terms: current W terms (for context, not used in GIO construction)

        Returns:
            list of monomial strings with R < max_R, or None if validation fails
        """
        # Determine required order from min R
        gauge_r = [r_charges.get(f, 1.0) for f in self.all_gauge_fields]
        if gauge_r:
            min_R = min(gauge_r)
            if min_R <= 0:
                return None  # non-positive R-charge, theory is inconsistent
            required_order = math.ceil(max_R / min_R)
        else:
            required_order = 4
        required_order = max(required_order, 4)  # minimum order 4

        # Build/extend if needed
        if required_order > self.max_order_built:
            if not self.build_and_validate(required_order):
                return None
            # Auto-save after building/extending
            self.save(self.cache_path())

        # Get products at the required order
        products = self.products_by_order.get(required_order,
                    self.products_by_order.get(self.max_order_built, []))

        # Filter by R < max_R
        def mono_R(mono):
            return sum(r_charges.get(f, 0) for f in mono.split('*'))

        gauge_ops = [p for p in products if mono_R(p) < max_R]

        # Add singlet fields and products
        all_ops = list(gauge_ops)
        for sf in self.singlet_fields:
            R_sf = r_charges.get(sf, 0)
            if R_sf > 0 and R_sf < max_R:
                all_ops.append(sf)
            # sf^2
            if R_sf > 0 and 2 * R_sf < max_R:
                all_ops.append(f"{sf}*{sf}")

        # Singlet × gauge ops
        for sf in self.singlet_fields:
            R_sf = r_charges.get(sf, 0)
            for gop in gauge_ops:
                if R_sf + mono_R(gop) < max_R:
                    all_ops.append(f"{sf}*{gop}")

        # Singlet × singlet (different)
        for i, sf1 in enumerate(self.singlet_fields):
            for sf2 in self.singlet_fields[i+1:]:
                if r_charges.get(sf1, 0) + r_charges.get(sf2, 0) < max_R:
                    all_ops.append(f"{sf1}*{sf2}")

        # Deduplicate by sorted field tuple
        seen = set()
        unique_ops = []
        for op in all_ops:
            key = tuple(sorted(op.split('*')))
            if key not in seen:
                seen.add(key)
                unique_ops.append(op)

        return unique_ops

    def save(self, path):
        """Save full GIO cache to JSON file."""
        data = {
            'lie_group': self.lie_group,
            'rep_names': self.rep_names,
            'rep_fields': self.rep_fields,
            'dynkin_labels': self.dynkin_labels,
            'singlet_fields': self.singlet_fields,
            'primitives': [
                {'degree': d, 'sym_type': s, 'rep': r, 'monomials': m}
                for d, s, r, m in (self.primitives or [])
            ],
            'products_by_order': {str(k): v for k, v in self.products_by_order.items()},
            'pe_counts_by_order': {
                str(k): {str(deg): cnt for deg, cnt in v.items()}
                for k, v in self.pe_counts_by_order.items()
            },
            'matching_by_order': self.matching_by_order,
            'max_order_built': self.max_order_built,
        }
        os.makedirs(os.path.dirname(path) if os.path.dirname(path) else '.', exist_ok=True)
        with open(path, 'w') as f:
            json.dump(data, f, indent=2)
        print(f"  GIO cache saved to {path}")

    @classmethod
    def load(cls, path):
        """Load GIO cache from JSON file."""
        with open(path) as f:
            data = json.load(f)

        cache = cls.__new__(cls)
        cache.lie_group = data['lie_group']
        cache.rep_names = data['rep_names']
        cache.rep_fields = data['rep_fields']
        cache.dynkin_labels = data['dynkin_labels']
        cache.singlet_fields = data['singlet_fields']

        cache.all_gauge_fields = []
        for name in cache.rep_names:
            cache.all_gauge_fields.extend(cache.rep_fields[name])

        cache.pe_labels = []
        cache.pe_mults = []
        cache._label_to_idx = {}
        for i, name in enumerate(cache.rep_names):
            dl = tuple(cache.dynkin_labels[i])
            if dl not in cache._label_to_idx:
                cache._label_to_idx[dl] = len(cache.pe_labels)
                cache.pe_labels.append(list(dl))
                cache.pe_mults.append(0)
            cache.pe_mults[cache._label_to_idx[dl]] += len(cache.rep_fields[name])

        cache.primitives = [
            (p['degree'], p['sym_type'], p['rep'], p['monomials'])
            for p in data.get('primitives', [])
        ]
        cache.products_by_order = {
            int(k): v for k, v in data.get('products_by_order', {}).items()
        }
        cache.pe_counts_by_order = {
            int(k): {eval(deg): cnt for deg, cnt in v.items()}
            for k, v in data.get('pe_counts_by_order', {}).items()
        }
        cache.matching_by_order = data.get('matching_by_order', {})
        cache.max_order_built = data.get('max_order_built', 0)

        print(f"  GIO cache loaded from {path} (max order {cache.max_order_built})")
        return cache

    def cache_path(self, base_dir='results'):
        """Standard path for this cache file."""
        matter_key = '_'.join(f"{n}{len(self.rep_fields[n])}" for n in self.rep_names)
        singlet_key = f"_s{len(self.singlet_fields)}" if self.singlet_fields else ""
        return os.path.join(base_dir, self.lie_group,
                           f"gio_{matter_key}{singlet_key}.json")


if __name__ == '__main__':
    import multiprocessing
    multiprocessing.set_start_method('fork')

    # Test: SU(2) 8 fund
    print("=== SU(2) 8 fund ===")
    cache = GIOCache(
        'A1',
        {'r0': [f'r0{i}' for i in range(1, 9)]},
        [[1]],
    )
    ok = cache.build_and_validate(order=4)
    print(f"  Validation: {'PASSED' if ok else 'FAILED'}")
    print(f"  Summary: {json.dumps(cache.summary(), indent=2)[:500]}")

    # Get operators with specific R-charges
    r = {f'r0{i}': 0.5 for i in range(1, 9)}
    ops = cache.get_operators(r, max_R=2.0)
    print(f"  Operators (R<2): {len(ops)}")

    # With lower R-charges (need higher order?)
    r2 = {f'r0{i}': 0.3 for i in range(1, 9)}
    ops2 = cache.get_operators(r2, max_R=2.0)
    print(f"  Operators (R<2, min_R=0.3): {len(ops2)}")

    # Test: SO(7) 3 adj
    print("\n=== SO(7) 3 adj ===")
    cache2 = GIOCache(
        'B3',
        {'phi': ['phi1', 'phi2', 'phi3']},
        [[0, 1, 0]],
    )
    ok2 = cache2.build_and_validate(order=4)
    print(f"  Validation: {'PASSED' if ok2 else 'FAILED'}")
    r3 = {'phi1': 2/3, 'phi2': 2/3, 'phi3': 2/3}
    ops3 = cache2.get_operators(r3, max_R=2.01)
    print(f"  Operators (R<2.01): {len(ops3)}")
