# Task 0.4: Documentation of Nf=2N.nb

## Overview

Mathematica notebook implementing the full classification pipeline for $\mathrm{SU}(N_c)$ gauge theories. Processes seed theories through iterative deformation up to depth 5, storing results in SQLite.

## Notebook Structure

| Section | Content |
|---------|---------|
| SU(N) central charges > a-maximization | `FindCharges` function |
| Superpotential to Graph | `toGlobalGraph` equivalence checker |
| Workers / DB | `workerFunction`, `SafeUpsert`, database setup |
| SU2nf4, SU2nf5, SU3nf5, SU3nf6, ... | Per-theory chapters, each with: Import data, Gauge Invariant Operators, Length0 through Length5, Distribution |

## Matter Field Convention

9 representation slots for $\mathrm{SU}(N_c)$:

| Index | Name | Dynkin index $\mu$ | Dimension $|R|$ |
|-------|------|-------------------|-----------------|
| 1 | X (singlet) | 0 | 1 |
| 2 | M (singlet) | 0 | 1 |
| 3 | q (fund) | 1/2 | $N_c$ |
| 4 | qb (anti-fund) | 1/2 | $N_c$ |
| 5 | phi (adjoint) | $N_c$ | $N_c^2 - 1$ |
| 6 | S (symmetric) | $(N_c+2)/2$ | $N_c(N_c+1)/2$ |
| 7 | Sb (anti-symmetric-bar) | $(N_c+2)/2$ | $N_c(N_c+1)/2$ |
| 8 | A (antisymmetric) | $(N_c-2)/2$ | $N_c(N_c-1)/2$ |
| 9 | Ab (antisymmetric-bar) | $(N_c-2)/2$ | $N_c(N_c-1)/2$ |

Multiplicity vector: $n = (N_X, N_M, N_q, N_{qb}, N_\phi, N_S, N_{Sb}, N_A, N_{Ab})$

Field variables are named by concatenating the name with an index: `q1`, `q2`, `M1`, `X1`, etc.

## FindCharges[nc, n, w]

Main function. Input: $N_c$, multiplicity vector $n$, superpotential term list $w$ (list of strings like `{"q1*q2", "M1*q1*q2"}`).

### Step 1: Setup
- Assigns symbolic R-charge variables `q[field]` to each field
- Defines Dynkin indices $\mu$ and dimensions for each representation slot

### Step 2: ABJ Anomaly
- Computes $\mu_{\mathrm{adj}} + \sum_i (q[\Phi_i] - 1) \mu_i = 0$
- Solves for one variable (the first appearing in the sum) and substitutes

### Step 3: Superpotential Constraints
- For each term in $w$: computes the total R-charge of the monomial
- If not already equal to 2, solves $R[W_k] = 2$ for one variable
- **Priority:** solves for $X$ or $M$ fields first (singlets), then others
- Each constraint eliminates one R-charge variable

### Step 4: Symmetry Identification
- Groups fundamental fields (`q1`, `q2`, ...) by their equivalence under `toGlobalGraph`
- Fields in the same equivalence class are assigned equal R-charges
- This further reduces the number of free variables

### Step 5: a-Maximization
- Constructs $\mathrm{Tr}\,R^3$ and $\mathrm{Tr}\,R$ from the R-charge variables
- Computes $a = (3/32)(3\,\mathrm{Tr}\,R^3 - \mathrm{Tr}\,R)$
- Checks: number of free variables = Total fields $-$ number of $W$ terms $-$ 1 (from ABJ)
- **Primary method:** `Solve[D[a, {vars}] == 0, vars]` (exact symbolic), with 1-hour timeout
- **Fallback:** `NSolve` with `WorkingPrecision -> 50`, also with 1-hour timeout
- If both fail: throws `"a-maximization timeout"`

### Step 6: Select Local Maximum
- Computes Hessian $H_{ij} = \partial^2 a / \partial x_i \partial x_j$
- Filters solutions where Hessian is negative definite (`NegativeDefiniteMatrixQ`)
- Selects the solution with the largest $a$ value (`bestSol`)
- If no local maximum exists: throws `"a has no local maximum"`

### Step 7: Consistency Checks
- **Real:** $a, c$ and all R-charges must be real numbers
- **Positive:** all R-charges $> 0$ and $\neq 2$; $a, c > 0$
- **Hofman-Maldacena:** $1/2 < a/c < 3/2$
- **Operator decoupling:** checks if any GIO has $R < 2/3$ (marks as `"operator decoupled"`)

### Step 8: Gauge-Invariant Operators
- **F-term relations:** computes $\partial W / \partial \Phi_i$ for all fields, then `GroebnerBasis[fterms, flatMat]`
- **GIO candidates:** combines flip fields ($M_j$, $X_j$) with `freeGIO` (the pre-computed list of gauge-invariant composites from matter fields)
- Also includes products: flip * GIO, flip * flip
- **Reduction:** `PolynomialReduce[gio, ringRelations, flatMat]` — keeps only the remainder (operators that are nonzero modulo F-terms and are monomials, not sums)
- **R-charge assignment:** `opcharges[expr]` computes $R[\mathcal{O}]$ by summing R-charges of constituent fields
- **Classification:**
  - `inequiv_relevant`: GIOs with $R < 2$, deduplicated by `toGlobalGraph`
  - `inequiv_flip`: GIOs with $R < 4/3$, deduplicated by `toGlobalGraph`
  - `all_relevant`, `all_flip`: full lists before deduplication

### Step 9: Output
Returns an Association with keys:
`theory`, `n`, `w`, `consistency`, `method`, `a`, `c`, R-charges per field, `rational` (exact form), `global` (flavor charges), `inequiv_relevant`, `inequiv_flip`, `all_relevant`, `all_flip`, `computation_time`

## toGlobalGraph[termList, maxM]

Converts a superpotential + matter content into a bipartite graph for equivalence checking.

- **Term nodes:** one per superpotential term (labeled `"1_Term"`, `"2_Term"`, ...)
- **Variable nodes:** one per field appearing in the terms
- **Edges:** connect term nodes to variable nodes, with multiplicity from powers
- **Chain structure:** handles repeated variables by creating chain nodes
- Used in `DeleteDuplicatesBy[..., toGlobalGraph]` to identify equivalent theories before running `FindCharges`

## freeGIO (Gauge-Invariant Operators)

**Hard-coded per theory.** For example, for SU(2) with 8 fundamentals (`SU2_8q`):
- `matter = {q1, q2, ..., q8}`
- `messon = Subsets[matter, {2}]` — all pairs $q_i q_j$ (mesons for SU(2))
- `quartic = Tuples[messon, 2]` — products of two mesons
- `freeGIO = Join[messon, quartic]`

This is the part that must be replaced by the Hilbert series for generalization.

## workerFunction[item]

Wrapper that calls `FindCharges[nc, item[[1]], item[[2]]]` and writes the result to the database via `SafeUpsert`.

## Iterative Pipeline (per theory chapter)

### Depth 0 (Length0)
```
length0 = Map[workerFunction, {{n_seed, {}}}]
```
Seed theory with $W = 0$.

### Depth n -> n+1 (Length n+1)
```
nw = for each consistent theory ii in length_n:
  - for each op in ii["inequiv_relevant"]:
      {ii["n"], Append[ii["w"], op]}                    # direct deformation
  - for each op in ii["inequiv_flip"]:
      {ii["n"] + {0,1,0,...}, Append[ii["w"], "M"<>...<>"*"<>op]}  # flip
```
Then:
1. Deduplicate `nw` via `ParallelMap[toGlobalGraph, ...]` + `DeleteDuplicatesBy`
2. Run `Map[workerFunction, rnw]` on unique candidates
3. Filter consistent results, accumulate into `total`
4. Export to `.wxf` file

### Accumulation
```
total = Join[length0, length1, ..., length5]
```

## Database (SQLite)

- Connected via `OpenSQLConnection[JDBC["SQLite", dbPath]]`
- `SafeUpsert[res]`: inserts or updates a row keyed by theory data
- Each theory chapter has its own database path

## Key Limitations for Generalization

1. **freeGIO is hard-coded** per gauge group and matter content (mesons, quartic products). Must be replaced by Hilbert series.
2. **Representation data ($\mu$, dim) is hard-coded** for $\mathrm{SU}(N_c)$ with the 9 fixed slots. Must be parameterized by Dynkin labels.
3. **`group` and `matter` functions in hilbert.py** are hard-coded for type $C$. Must accept arbitrary Lie types.
4. **Parallelization:** `workerFunction` uses `ParallelMap` for deduplication but sequential `Map` for `FindCharges`. Could benefit from parallelization at the `FindCharges` level too.
