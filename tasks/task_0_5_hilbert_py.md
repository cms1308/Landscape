# Task 0.5: Documentation of hilbert.py

## Overview

Python script that computes the Hilbert series for a gauge theory using:
1. **LiE** — computes the plethystic exponential (gauge singlet projection) via symmetric tensor products and the Frobenius formula
2. **Mathematica** (`wolframscript`) — computes the plethystic logarithm (Möbius inversion) to extract chiral ring generators and relations
3. **multiprocessing** — parallelizes the PE computation over partitions

## Dependencies

- `lie` (LiE computer algebra system for Lie groups) — called via `subprocess`
- `wolframscript` (Mathematica kernel) — called via `subprocess.Popen`
- Python standard library: `multiprocessing`, `subprocess`, `ast`, `itertools`

## Input Specification

Defined in `__main__`:

```python
rank = 4                                    # Lie algebra rank
count = [0, 0, 0, 0, 0, 1, 0, 1, 0]       # [X, M, q, qb, adj, s, sb, a, ab]
order = 6                                   # truncation order in fugacities
```

- `count`: multiplicity vector using the same 9-slot convention as `Nf=2N.nb`
- `order`: maximum total degree in the fugacity expansion

## Functions

### `group(rank)` — Gauge group specification

```python
def group(rank):
    return 'C' + str(rank)
```

**Hard-coded to type C (Sp).** Returns LiE group name, e.g. `"C4"` for $\mathrm{Sp}(4)$.

**Must generalize:** accept Lie type as parameter, return `"A4"`, `"B3"`, `"D5"`, `"G2"`, `"F4"`, `"E6"`, etc.

### `matter(rank, num)` — Representations with per-copy fugacities

Builds a list of Dynkin label vectors, one per copy of each representation. For `count = [0,0,2,0,1,0,0,0,0]` (2 fundamentals + 1 adjoint of $C_r$), returns:

```
[[1,0,...,0], [1,0,...,0], [1,0,...,0,1]]
```

Each copy gets its own entry → each gets its own fugacity $z_{i,a}$ in the PE.

**Hard-coded Dynkin labels for type C:**

| Slot | Rank $\geq 2$ | Rank 1 |
|------|---------------|--------|
| X, M (singlet) | $[0,\ldots,0]$ | $[0]$ |
| fund | $[1,0,\ldots,0]$ | $[1]$ |
| antifund | $[0,\ldots,0,1]$ | $[1]$ |
| adj | $[1,0,\ldots,0,1]$ | $[2]$ |
| sym | $[2,0,\ldots,0]$ | $[2]$ |
| symbar | $[0,\ldots,0,2]$ | $[2]$ |
| anti | $[0,1,0,\ldots,0]$ | $[0]$ |
| antibar | $[0,\ldots,0,1,0]$ | $[0]$ |

**Must generalize:** accept arbitrary Dynkin labels, not restricted to 9 fixed slots or type C.

### `matters2(rank, num)` — Representations with per-type fugacities

Same as `matter` but includes each representation type **once** (not per copy). Used with `num2` which passes the multiplicities separately. Each representation type gets one fugacity $z_i$.

### `num(nums)` / `num2(nums)` — Multiplicity formatting

- `num(nums)`: returns `[1,1,1,...,1]` of length `sum(nums)` — each copy counted separately
- `num2(nums)`: returns only nonzero entries from `nums` — multiplicities per representation type

These correspond to the two fugacity modes (per-copy vs per-type).

### `PE(count, matters, parts, nums, rank)` — Plethystic exponential via LiE

Core computation. Generates a LiE script and runs it via `subprocess`.

**LiE algorithm** (for a single partition vector `part[j]`):

```
For each representation k:
  If part[j,k] == 0:
    tensor with identity (no contribution)
  Else:
    Compute Frobenius character: frob = from_part(partitions(part[j,k]))
    Filter by multiplicity: keep only terms where sum <= num[k]
    For each valid Frobenius term:
      Build the representation by:
        rep = product over columns n of:
          C(count, ffrob[m,n]) * p_tensor(ffrob[m,n], sym_tensor(n, matter[k], group), group)
    Sum all such reps
  Tensor all representations together
If the result is a singlet (coef of trivial rep != 0):
  Record (partition, coefficient)
```

Key LiE functions used:
- `sym_tensor(n, R, G)`: $n$-th symmetric power of representation $R$ under group $G$
- `p_tensor(n, R, G)`: $n$-th Adams (plethystic) power
- `tensor(R1, R2, G)`: tensor product projected onto $G$-representations
- `from_part(partitions(n))`: Frobenius character from integer partitions
- `poly_one(d)`: trivial representation (identity in the polynomial ring)
- `coef(term, 1)`: coefficient of the trivial representation

**Output:** list of `[partition_vector, multiplicity]` pairs — the gauge-singlet content at each order.

**Parallelization:** each partition is processed independently. `run`/`run2` dispatch all partitions via `mp.Pool().starmap(PE, ...)`.

### `PL(strs)` — Plethystic logarithm (per-copy fugacities)

Takes the PE output and runs Mathematica to compute:

$$\mathrm{PL}[f(z)] = \sum_{k=1}^{K} \frac{\mu(k)}{k} \log f(z^k)$$

where $\mu(k)$ is the Möbius function and the substitution $z \to z^k$ is applied to all fugacities.

**Fugacity assignment:** field $a$ of representation type $i$ gets fugacity $z_{i,a}$ (subscripted by both type and copy index). This distinguishes individual copies.

The result is expanded as a power series in all fugacities up to `order`, then `Normal // Expand`.

**Output:** two lines printed — the PE as a polynomial, then the PL.
- Positive terms in PL = **generators** of the chiral ring
- Negative terms in PL = **relations** among generators

### `PL2(strs)` — Plethystic logarithm (per-type fugacities)

Same as `PL` but with one fugacity $z_i$ per representation type (not per copy). Used with `matters2`/`num2`.

### `run(matters, parts, nums, rank)` / `run2(...)`

End-to-end pipeline:
1. Parallel PE over all partitions via `mp.Pool().starmap`
2. Filter out `None` results (partitions with no singlet content)
3. Format PE output as Mathematica-compatible string
4. Call `PL` or `PL2` to compute the plethystic logarithm

`run` uses per-copy fugacities (`matter` + `num` + `PL`).
`run2` uses per-type fugacities (`matters2` + `num2` + `PL2`).

### `Hilbert(strs)` — (Unused/broken)

Appears to be an incomplete function — has a typo (`proc.communite()`) and a stray `ng1` in the Mathematica code. Not called anywhere.

## Example Run

```python
rank = 4
count = [0, 0, 0, 0, 0, 1, 0, 1, 0]  # 1 symmetric + 1 antisymmetric of Sp(4)
order = 6
```

This computes the Hilbert series for $\mathrm{Sp}(4)$ with 1 symmetric $[2,0,0,0]$ and 1 antisymmetric $[0,1,0,0]$, expanded to order 6. Output written to `Sp4s1a1nf0.txt`.

## Commented-Out Example

```python
rank = 6
count = [0, 0, 0, 0, 0, 2, 2, 0, 0]  # 2 symmetric + 2 anti-symmetric
```

Uses `run` (per-copy fugacities) instead of `run2` — output to `SU5s2S2nf0_1.txt`. Note the filename says `SU5` but the code uses type C (Sp).

## What the PE Output Means

Each entry `[n1, n2, ..., n_k, coeff]` means: the gauge-singlet sector contains `coeff` copies of operators at fugacity degree $(n_1, n_2, \ldots, n_k)$.

After PL:
- Positive coefficient at degree $(n_1, \ldots, n_k)$: there are that many **independent generators** at that degree
- Negative coefficient: **relations** (syzygies) at that degree

For our project, we need the full PE (not just PL) to enumerate all GIOs at each R-charge level.

## Limitations for Generalization

1. `group(rank)` returns only type $C$
2. `matter`/`matters2` have 9 hard-coded representation slots with type-$C$ Dynkin labels
3. No F-term reduction — the Hilbert series is computed for the free theory (no superpotential). F-terms must be imposed separately via `GroebnerBasis` as in `Nf=2N.nb`
4. `Hilbert` function is broken (typo + stray code)
5. Notification uses Slack webhook (redacted); could be updated to use Telegram or removed
