"""
Generalized Hilbert series computation for arbitrary simple Lie groups.

Computes the plethystic exponential (PE) of gauge-invariant operators
using LiE for group theory and multiprocessing for parallelism.
Optionally computes the plethystic logarithm (PL) via Mathematica.

Generalization of references/hilbert.py:
  - group: accepts any Lie type (A, B, C, D, G2, F4, E6, E7, E8)
  - matter: accepts arbitrary Dynkin labels, not restricted to 9 fixed slots
  - multiplicities: arbitrary
"""

import subprocess
import multiprocessing as mp
import ast
import time
import sys
import os
import json
from itertools import product, repeat


def PE(count, lie_group, matters_str, nums_str, parts):
    """Compute plethystic exponential for a single partition via LiE.

    Args:
        count: index into parts list
        lie_group: LiE group name (e.g. "A2", "B3", "G2")
        matters_str: LiE-formatted matter representation list
        nums_str: LiE-formatted multiplicity list
        parts: list of all partitions
    Returns:
        String of LiE output (singlet content) or None
    """
    lcode = """
    maxnodes 9999999
    maxobjects 9999999
    group=%s
    matter=%s
    num=%s
    part=%s
    comb(int n,m)={
    class_ord([n+1])/(class_ord([n-m+1])*class_ord([m+1]))
    }
    single=poly_one(n_cols(matter))
    sol=null(0,n_cols(part)+1)
        for j=1 to n_rows(part) do
            term=single;
            for k=1 to size(part[j]) do
                if part[j,k]==0 then term=tensor(term,single,group);
                else frob=from_part(partitions(part[j,k]));
                    ffrob=null(0,n_cols(frob)+1);
                    for l=1 to n_rows(frob) do
                        if all_one(part[j,k])*(frob[l]+(l/(n_rows(frob))))<=num[k] then ffrob=ffrob+(frob[l]+(l/(n_rows(frob)))); fi;
                    od;
                zero=0*single;
                for m=1 to n_rows(ffrob) do
                    rep=single;
                    count=num[k];
                    for n=1 to n_cols(ffrob) do
                        rep=comb(count,ffrob[m,n])*tensor(rep,p_tensor(ffrob[m,n], sym_tensor(n,matter[k],group),group),group);
                        count=count-ffrob[m,n];
                    od;
                    zero=zero+rep;
                od;
                term=tensor(term,zero,group);
                fi;
            od;
            if coef(term,1)!=0 && term[1]/coef(term,1)==single then sol=sol+(part[j]+coef(term,1)); fi;
        od;
    print(sol);
    """ % (lie_group, matters_str, nums_str, '[' + str(parts[count]) + ']')

    res = subprocess.run(
        'lie', shell=True, input=lcode,
        capture_output=True, encoding='UTF-8'
    ).stdout[53:].replace('\n', '').replace(" ", "")

    if 'null' in res:
        return None
    else:
        return res


def PL_mathematica(pe_str, counts_str, order):
    """Compute plethystic logarithm via Mathematica (per-type fugacities).

    Args:
        pe_str: Mathematica-formatted PE output
        counts_str: Mathematica-formatted multiplicity list {n1, n2, ...}
        order: truncation order
    Returns:
        String with PE and PL as Mathematica expressions
    """
    mcode = """
    count = %s; count2 = {}; count3 = {};count4 = 1;
    For[i = 1, i <= Length[count], i++,
        If[count[[i]]!=0, count2=Append[count2, count4]; count4++;];
    ];

    seriesvar = Subscript[z, count2[[#]]] & /@ Range[Length[count2]];
    f = Total[#[[Length[#]]] Product[Subscript[z, count2[[i]]]^#[[i]], {i, 1,Length[#] - 1}] & /@ %s];
    Print[f];
    PL[x_, k_] := Sum[(MoebiusMu[i] Log[x])/i /.Subscript[z,w_] :>  Subscript[z,w]^i, {i, 1, k}];
    pl = PL[f, %d];
    For[i = 1, i <= Length[%s[[1]]]-1, i++,pl = Series[pl, {seriesvar[[i]], 0, %d}] // Normal // Expand;];
    Print[pl];
    """ % (counts_str, pe_str, order, pe_str, order)

    proc = subprocess.Popen(
        ['wolframscript', '-code', mcode],
        stdout=subprocess.PIPE
    )
    (out, err) = proc.communicate()
    return out.decode('ascii').replace("Null", "")


def compute_hilbert_series(lie_group, dynkin_labels, multiplicities, order,
                           compute_pl=True):
    """Compute Hilbert series for a gauge theory.

    Args:
        lie_group: LiE group name (e.g. "A2", "B3", "G2", "F4")
        dynkin_labels: list of Dynkin label lists, one per rep type
                       e.g. [[1,0], [0,1], [1,1]] for fund, antifund, adj of A2
                       Singlets should be omitted (they don't participate in gauge projection)
        multiplicities: list of integers, number of copies of each rep type
        order: truncation order for the series expansion
        compute_pl: if True, also compute the plethystic logarithm

    Returns:
        dict with keys:
            'pe_raw': list of (degree_vector, multiplicity) pairs
            'pl': string of PL expression (if compute_pl)
            'time': computation time in seconds
    """
    start = time.time()

    # Filter out singlets (Dynkin label all zeros) — they don't enter gauge projection
    gauge_labels = []
    gauge_mults = []
    for dl, m in zip(dynkin_labels, multiplicities):
        if any(x != 0 for x in dl):  # non-singlet
            gauge_labels.append(dl)
            gauge_mults.append(m)

    if not gauge_labels:
        # All singlets — no gauge projection needed, trivial Hilbert series
        return {'pe_raw': [], 'pl': '1', 'time': 0}

    # Format for LiE
    matters_str = str(gauge_labels).replace(" ", "")
    nums_str = str(gauge_mults).replace(" ", "")
    counts_str = str(gauge_mults).replace("[", "{").replace("]", "}")

    # Generate all partitions up to given order
    n_types = len(gauge_mults)
    parts = list(product(range(0, order + 1), repeat=n_types))
    # Convert to LiE format
    parts_list = [list(p) for p in parts]

    print(f"  Computing PE: {lie_group}, {len(gauge_labels)} rep types, "
          f"{len(parts_list)} partitions, order={order}")

    # Parallel PE computation
    with mp.Pool() as pool:
        ans = pool.starmap(
            func=PE,
            iterable=zip(
                range(len(parts_list)),
                repeat(lie_group),
                repeat(matters_str),
                repeat(nums_str),
                repeat(parts_list)
            )
        )

    # Filter None results
    pe_results = [r for r in ans if r is not None]

    # Format for Mathematica
    pe_str = str(pe_results).replace(
        "[[", "{").replace("]]", "}").replace(
        "[", "{").replace("]", "}").replace("'", "").replace(" ", "")

    elapsed_pe = time.time() - start
    print(f"  PE done: {len(pe_results)} non-zero terms ({elapsed_pe:.1f}s)")

    result = {
        'pe_raw': pe_results,
        'time_pe': elapsed_pe,
    }

    # Compute PL
    if compute_pl and pe_results:
        print(f"  Computing PL via Mathematica...")
        pl_output = PL_mathematica(pe_str, counts_str, order)
        result['pl'] = pl_output
        result['time_total'] = time.time() - start
        print(f"  PL done ({result['time_total']:.1f}s total)")
    else:
        result['pl'] = None
        result['time_total'] = elapsed_pe

    return result


if __name__ == '__main__':
    # Example: G2 with 1 fundamental [1,0] + 1 adjoint [0,1]
    if len(sys.argv) > 1:
        # Usage: python hilbert.py <config.json>
        with open(sys.argv[1]) as f:
            config = json.load(f)
        lie_group = config['group']
        dynkin_labels = config['dynkin_labels']
        multiplicities = config['multiplicities']
        order = config.get('order', 6)
    else:
        # Default test case
        lie_group = "G2"
        dynkin_labels = [[1, 0], [0, 1]]
        multiplicities = [1, 1]
        order = 6

    print(f"Hilbert series for {lie_group}")
    print(f"  Reps: {dynkin_labels}")
    print(f"  Mults: {multiplicities}")
    print(f"  Order: {order}")

    result = compute_hilbert_series(
        lie_group, dynkin_labels, multiplicities, order
    )

    print(f"\nPE terms: {len(result['pe_raw'])}")
    if result['pl']:
        print(f"\nPL output:\n{result['pl'][:500]}")
