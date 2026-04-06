# Task 1.3: Gauge-Invariant Operators from the Hilbert Series

## Problem

In `Nf=2N.nb`, the gauge-invariant operators (GIOs) are constructed by hand for each gauge group and matter content. For example, for SU(2) with 8 fundamentals:

```mathematica
matter = {q1, q2, ..., q8};
messon = Subsets[matter, {2}];          (* pairs q_i q_j *)
quartic = Tuples[messon, 2];            (* products of two mesons *)
freeGIO = Join[messon, quartic];
```

This is group-specific and cannot be generalized to arbitrary gauge groups and representations.

## Solution: Plethystic Exponential (PE)

The Hilbert series, computed via the PE in `hilbert.py`, automatically enumerates all gauge-invariant operators for **any** gauge group and matter content. No hand-coded invariant tensor analysis is needed.

### What the PE computes

For gauge group $G$ with matter fields $\Phi_i^{(a)}$ ($a = 1, \ldots, N_i$) in representations $R_i$, the PE of the single-particle partition function gives the **multi-particle gauge-invariant partition function**:

$$H(z_1, z_2, \ldots) = \mathrm{PE}\left[\sum_i \sum_{a=1}^{N_i} z_{i,a}\, \chi_{R_i}\right]_G$$

where $z_{i,a}$ are fugacities for each field, $\chi_{R_i}$ is the character of $R_i$, and the subscript $G$ means projection onto gauge singlets.

Expanding $H$ in powers of the fugacities:

$$H = 1 + \sum_{\vec{n}} c(\vec{n})\, z_1^{n_1} z_2^{n_2} \cdots$$

Each term $c(\vec{n})\, z_1^{n_1} z_2^{n_2} \cdots$ with $c(\vec{n}) \neq 0$ represents gauge-invariant operators built from $n_1$ copies of field 1, $n_2$ copies of field 2, etc. The coefficient $c(\vec{n})$ counts the number of independent such operators.

### What we extract

At each fixed point, after $a$-maximization determines the R-charges $r_i$ of each field, we substitute $z_{i,a} \to t^{r_i}$ (or keep individual fugacities for flavor resolution) and read off:

- **All GIOs with $R < 2$** (relevant): candidates for direct deformation $\delta W = \mathcal{O}$
- **All GIOs with $R < 4/3$** (super-relevant): candidates for flip $\delta W = M \cdot \mathcal{O}$
- **All GIOs with $R < 2/3$** (unitarity-violating): trigger decoupling

The truncation order in `hilbert.py` controls how high in R-charge we enumerate.

### PE vs PL

| | PE (Hilbert series) | PL (Plethystic logarithm) |
|---|---|---|
| Gives | All GIOs at each degree | Generators + relations of chiral ring |
| Includes composites | Yes | No (only primitive generators) |
| What we need | **This one** — to find all relevant operators | Useful as cross-check |

A composite of two generators, each with $R < 1$, is itself a relevant operator with $R < 2$. The PL would miss it since it only lists generators. Therefore we use the **PE** for operator enumeration.

## Integration with F-term Reduction

The PE gives the GIOs of the **free theory** (no superpotential). When $W \neq 0$, some GIOs become zero in the chiral ring due to F-term relations $\partial W / \partial \Phi_i = 0$.

The pipeline:
1. **PE** (via `hilbert.py` / LiE): enumerate all gauge-invariant operators up to the relevant R-charge
2. **F-term reduction** (via Mathematica `GroebnerBasis` / `PolynomialReduce`, as in `Nf=2N.nb`): impose $\partial W / \partial \Phi_i = 0$ and remove operators that are zero in the chiral ring
3. **R-charge assignment**: assign $R[\mathcal{O}] = \sum (\text{constituent R-charges})$ from $a$-maximization
4. **Classification**: sort into relevant ($R < 2$), super-relevant ($R < 4/3$), unitarity-violating ($R < 2/3$)
5. **Deduplication**: use `toGlobalGraph` to identify equivalent operators under flavor permutations

## How This Replaces freeGIO

| `Nf=2N.nb` (old) | Hilbert series (new) |
|---|---|
| Hand-code `matter`, `messon`, `quartic` | PE automatically computes all singlets |
| Group-specific (SU(2) mesons = $q_i q_j$) | Group-independent (works for any $G$, any reps) |
| Fixed operator degree (mesons = degree 2, quartic = degree 4) | Systematic to any degree (controlled by `order` parameter) |
| May miss operators (e.g. baryons, dressed mesons) | Complete within the truncation order |

## Practical Workflow

For a theory with gauge group $G$ (LiE type + rank), matter content specified by Dynkin labels and multiplicities:

```
1. Call hilbert.py with:
   - group = LiE type (e.g. "A2", "B3", "G2")
   - representations = list of Dynkin labels
   - multiplicities = number of copies of each
   - order = truncation (set high enough to capture all R < 2 operators)

2. Parse PE output: list of (degree_vector, multiplicity) pairs

3. Map degree vectors to symbolic operators (products of field variables)

4. If W != 0: reduce modulo F-terms via GroebnerBasis

5. Assign R-charges, classify, deduplicate
```

## Truncation Order

The PE is computed as a power series truncated at `order`. To capture all operators with $R[\mathcal{O}] < 2$, the order must satisfy:

$$\texttt{order} \geq \left\lfloor \frac{2}{\min_i r_i} \right\rfloor$$

where $r_i$ are the R-charges of the elementary fields. For seed theories ($W = 0$), the R-charges are determined by $a$-maximization and are typically $O(1)$, so `order = 6` is usually sufficient. For deeply deformed theories, R-charges of some fields may be smaller, requiring higher order.
