"""
Filter asymptotically free matter content for interacting IR fixed points.
Runs FindCharges (a-maximization) on each candidate from seed_matters.json.
"""

import json
import subprocess
import sys
import time
from fractions import Fraction

AMAX_PATH = "src/amax/FindCharges.wl"

# Group data: {group: (dimG, Tadj)}
GROUP_DATA = {
    'A1': (3, 2), 'A2': (8, 3), 'A3': (15, 4), 'A4': (24, 5),
    'B2': (10, 3), 'B3': (21, 5),
    'C2': (10, 3), 'C3': (21, 4),
    'D4': (28, 6),
    'G2': (14, 4), 'F4': (52, 9),
}

# Map from Dynkin labels to (name_prefix, T_R, dim_R) for each group
# Load from all_reps.json
def load_rep_lookup(all_reps, group):
    """Build lookup: tuple(dynkin_label) -> {T_R, dim, reality, name}"""
    lookup = {}
    for r in all_reps[group]['reps']:
        label = tuple(r['dynkin_label'])
        lookup[label] = {
            'T_R': r['T_R'],
            'dim': r['dim'],
            'reality': r['reality'],
            'name': r.get('name', str(label)),
        }
        # Also add conjugate
        if r['contragr']:
            conj = tuple(r['contragr'])
            lookup[conj] = {
                'T_R': r['T_R'],  # same T
                'dim': r['dim'],  # same dim
                'reality': r['reality'],
                'name': r.get('name', str(label)) + 'b',
            }
    return lookup


def make_name(label, rep_idx, is_conj=False):
    """Generate a short field name from Dynkin label."""
    s = 'r' + str(rep_idx)
    if is_conj:
        s += 'b'
    return s


def build_repinfo(matter_entry, rep_lookup, group):
    """Convert a seed matter entry to FindCharges repInfo format.

    matter_entry: list of [dynkin_label, N_copies, reality, dim]
    Returns: list of {name, T_R, dim_R, N_copies} tuples for Mathematica
    """
    repinfo = []
    for idx, (label, n, reality, dim_r) in enumerate(matter_entry):
        label_t = tuple(label)
        T_R = rep_lookup[label_t]['T_R']

        if reality == 'complex':
            # R + R_bar: add both as separate rep types
            name_r = make_name(label, idx)
            name_rb = make_name(label, idx, is_conj=True)
            repinfo.append((name_r, T_R, dim_r, n))
            repinfo.append((name_rb, T_R, dim_r, n))  # conjugate has same T, dim
        else:
            name_r = make_name(label, idx)
            repinfo.append((name_r, T_R, dim_r, n))

    return repinfo


def run_findcharges(dimG, Tadj, repinfo, w=[]):
    """Call FindCharges via wolframscript and parse the result."""
    # Format repInfo for Mathematica
    ri_str = '{' + ', '.join(
        f'{{"{name}", {T_R}, {dim_r}, {n}}}'
        for name, T_R, dim_r, n in repinfo
    ) + '}'

    w_str = '{' + ', '.join(f'"{t}"' for t in w) + '}'

    mcode = f"""
Get["{AMAX_PATH}"];
result = FindCharges[{dimG}, {Tadj}, {ri_str}, {w_str}];
Print[result["consistency"]];
Print[N[result["a"], 15]];
Print[N[result["c"], 15]];
Do[
  Print[k, " -> ", N[result[k], 15]],
  {{k, Select[Keys[result],
    StringQ[#] && !MemberQ[{{"n","w","consistency","method","a","c",
      "rational","computation_time","global"}}, #] &]}}
];
"""
    proc = subprocess.run(
        ['wolframscript', '-code', mcode],
        capture_output=True, encoding='UTF-8', timeout=120
    )
    lines = [l.strip() for l in proc.stdout.split('\n') if l.strip()]

    if len(lines) < 3:
        return {'consistency': 'error: ' + proc.stderr[:200]}

    result = {
        'consistency': lines[0],
        'a': lines[1],
        'c': lines[2],
    }
    # Parse R-charges
    for line in lines[3:]:
        if ' -> ' in line:
            k, v = line.split(' -> ', 1)
            result[k.strip()] = v.strip()

    return result


def filter_group(group, seeds_data, all_reps):
    """Run a-maximization on all seed candidates for a group."""
    dimG, Tadj = GROUP_DATA[group]
    rep_lookup = load_rep_lookup(all_reps, group)
    seeds = seeds_data[group]['seeds']

    print(f"\n{'='*70}")
    print(f"  {group}  |  {len(seeds)} candidates")
    print(f"{'='*70}")

    consistent = []
    inconsistent_counts = {}

    for i, seed in enumerate(seeds):
        matter = seed['matter']
        repinfo = build_repinfo(matter, rep_lookup, group)
        desc = seed['description']

        result = run_findcharges(dimG, Tadj, repinfo)
        status = result['consistency']

        if status == 'consistent':
            consistent.append({
                'description': desc,
                'matter': matter,
                'a': result.get('a', ''),
                'c': result.get('c', ''),
                'rcharges': {k: v for k, v in result.items()
                             if k not in ('consistency', 'a', 'c')},
            })
            print(f"  [{i+1}/{len(seeds)}] CONSISTENT  a={result.get('a','?')[:8]}  {desc}")
        else:
            inconsistent_counts[status] = inconsistent_counts.get(status, 0) + 1
            if (i+1) % 20 == 0:
                print(f"  [{i+1}/{len(seeds)}] ... ({len(consistent)} consistent so far)")

    print(f"\n  Summary: {len(consistent)} consistent / {len(seeds)} total")
    for reason, count in sorted(inconsistent_counts.items()):
        print(f"    {reason}: {count}")

    return consistent


if __name__ == '__main__':
    with open('data/all_reps.json') as f:
        all_reps = json.load(f)
    with open('data/seed_matters.json') as f:
        seeds_data = json.load(f)

    groups = sys.argv[1:] if len(sys.argv) > 1 else sorted(GROUP_DATA.keys())

    all_results = {}
    for grp in groups:
        if grp not in seeds_data:
            continue
        start = time.time()
        consistent = filter_group(grp, seeds_data, all_reps)
        elapsed = time.time() - start
        print(f"  Time: {elapsed:.1f}s")
        all_results[grp] = consistent

    with open('data/consistent_seeds.json', 'w') as f:
        json.dump(all_results, f, indent=2)
    print(f"\nSaved to data/consistent_seeds.json")
