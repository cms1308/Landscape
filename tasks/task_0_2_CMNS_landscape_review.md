# Task 0.2: Review of Large Landscape of 4d SCFTs (2408.02953)

**Paper:** M. Cho, K. Maruyoshi, E. Nardoni, J. Song, "Large Landscape of 4d Superconformal Field Theories from Small Gauge Theories"

## Overview

Systematic classification of 4d N=1 SCFTs starting from rank-1 and rank-2 Lagrangian gauge theories. By iteratively applying relevant deformations and flippings, the authors identify **7,346 inequivalent fixed points** cataloged at `https://qft.kaist.ac.kr/landscape`.

## Seed Theories

Starting points: Lagrangian gauge theories at conformal fixed points with W = 0.

| Seed | Gauge Group | Matter |
|------|-------------|--------|
| 1 | SU(2) | 1 adjoint |
| 2 | SU(2) | 2 adjoints |
| 3 | SU(3) | 1 adjoint |
| 4 | SO(5) | 1 adjoint |
| 5 | SO(5) | 2 adjoints |
| 6 | SO(5) | 1 traceless symmetric + 1 vector |
| 7 | SO(5) | 1 traceless symmetric + 2 vectors |
| 8 | Sp(2) | 1 rep in 14 + 2 fundamentals |
| 9 | Sp(2) | 1 adjoint |
| 10 | G_2 | 1 adjoint |

These are chosen to be inside the conformal window. All have simple gauge groups of rank 1 or 2 with small matter content.

## The Deformation Algorithm

### Step 1: Operator Enumeration

At each fixed point, enumerate all gauge-invariant chiral operators (GIOs) and compute their R-charges using $a$-maximization.

Classify:
- **Relevant:** $R[\mathcal{O}] < 2$ (can appear as $\delta W = \mathcal{O}$)
- **Super-relevant:** $R[\mathcal{O}] < 4/3$ (can be flipped: $\delta W = M \cdot \mathcal{O}$)

### Step 2a: Direct Deformation

For each relevant operator $\mathcal{O}$:
- Add $\delta W = \mathcal{O}$ to the superpotential
- Impose $R[\mathcal{O}] = 2$ as a new constraint
- Redo $a$-maximization

**Special case — single U(1):** If $\mathcal{O}$ is charged under a single U(1) flavor symmetry, the constraint $R_{\text{IR}}(\mathcal{O}) = 2$ uniquely fixes the mixing parameter: $\epsilon = (2 - R(\mathcal{O})) / J_1(\mathcal{O})$, with no optimization needed.

### Step 2b: Flip Deformation

For each super-relevant operator $\mathcal{O}$ with $R[\mathcal{O}] < 4/3$:
- Introduce gauge-singlet chiral field $M$
- Add $\delta W = M \cdot \mathcal{O}$
- The flip field gets $R[M] = 2 - R[\mathcal{O}] > 2/3$ (healthy)
- The original operator $\mathcal{O}$ becomes irrelevant ($R > 2$) in the new theory

**Central charge shift from a pure flip:**

$$\delta a = \frac{3}{32}(1 - R[\mathcal{O}])(3R[\mathcal{O}]^2 - 6R[\mathcal{O}] + 2)$$

$$\delta c = \frac{1}{32}(1 - R[\mathcal{O}])(9R[\mathcal{O}]^2 - 18R[\mathcal{O}] + 4)$$

The slope in the $(a,c)$ plane: $s = (9R^2 - 18R + 4)/(9R^2 - 18R + 6)$, ranging from $\sim 5/3$ to 2 for $2/3 < R < 4/3$.

**Why $R < 4/3$?** The flip field $M$ has $R[M] = 2 - R[\mathcal{O}]$. For $M$ to satisfy the unitarity bound, we need $R[M] \geq 2/3$, i.e. $R[\mathcal{O}] \leq 4/3$. If $R[\mathcal{O}] = 4/3$ exactly, $R[M] = 2/3$ and $M$ is a free field (marginal case).

### Step 3: Operator Decoupling

After $a$-maximization, if any operator $\mathcal{O}_d$ has $R[\mathcal{O}_d] < 2/3$:
- It violates unitarity and must decouple as a free field
- Introduce flip field $X$ with $\delta W = X \cdot \mathcal{O}_d$, enforcing $R[\mathcal{O}_d] = 2/3$
- Subtract free-field contribution from central charges:

$$\delta a = \frac{3}{32}\left[\frac{2}{9} - \left(3(R-1)^3 - (R-1)\right)\right]$$

This counts as a superpotential deformation (increments depth).

### Step 4: Iterate

Apply Steps 1–3 to every new fixed point. Stop when no new non-trivial fixed points are generated. The iteration depth defines "generations."

## Consistency Checks (Section 2.2)

Every candidate fixed point is tested against:

1. **Unitarity:** $R[\mathcal{O}] \geq 2/3$ for all chiral primaries (after decoupling)
2. **Hofman-Maldacena:** $1/2 \leq a/c \leq 3/2$
3. **Positive central charges:** $a, c > 0$
4. **$a$-theorem:** $a_{\text{child}} < a_{\text{parent}}$

### Superconformal Index Tests

The paper uses the **superconformal index** for detailed checks:

$$\mathcal{I}(t,y;x) = \text{Tr}\,(-1)^F t^{3(R+2j_1)} y^{2j_2} x^f$$

**Reduced index** (strips derivative descendants):

$$\hat{\mathcal{I}} = (1-t^3 y)(1-t^3/y)(\mathcal{I}-1)$$

From the reduced index:
- **Relevant operators:** coefficient of $t^{3R}$ (with $R < 2$) at $y^0$
- **Marginal operators:** coefficient of $t^6$ minus conserved currents gives $\alpha$ (conformal manifold dimension)
- **Conserved currents:** coefficient of $t^6$ from adjoint of flavor symmetry

### N=2 Enhancement Detection

- **Sufficient:** positive coefficient of $t^7(y + 1/y)$
- **N=3:** coefficient 2 for $t^7(y+1/y)$

**Note for our project:** We use the Hilbert series instead of the superconformal index, so these index-level tests (N=2 enhancement, conserved currents, conformal manifold dimension) are NOT available to us.

## Equivalence Criteria

Two fixed points are **inequivalent** if they differ in:
- Central charges $(a,c)$
- Operator spectrum (sorted R-charges of GIOs)
- Superconformal index (when computed)

Theories reached by different paths can be identified if all invariants match.

**Non-commutativity:** Different orderings of deformations can lead to different fixed points — the landscape is not a simple lattice.

## Key Results

| Quantity | Value |
|----------|-------|
| Total fixed points | 7,346 |
| $a/c$ range | 0.7228 to 1.2100 |
| Minimal SCFT | $(a,c) = (633/2000, 683/2000)$ |
| Theories with $a = c$ | Found (holographic relevance) |
| Conformal manifolds ($\alpha > 0$) | Found |
| Possible new N=2 SCFT | 1 candidate |

## What Our Project Generalizes

| Aspect | CMNS | Our Project |
|--------|------|-------------|
| Gauge groups | 5 specific (SU(2), SU(3), SO(5), Sp(2), G_2) | All simple groups |
| Rank | 1–2 | Arbitrary |
| Matter | Small, specific reps | All reps with $T(R) \leq T_{\max}$ |
| GIO enumeration | Superconformal index | Hilbert series + F-terms |
| Emergent symmetry | Detectable via index | Not detectable |
| Seed theories | 10 hand-picked | Systematically enumerated |

## Notation Adopted from CMNS

- R-charges: $R[\mathcal{O}]$, with $\Delta = (3/2)R$
- Relevant: $R < 2$; super-relevant: $R < 4/3$
- Flip: $\delta W = M \cdot \mathcal{O}$, with $R[M] = 2 - R[\mathcal{O}]$
- Depth = number of superpotential terms
- Central charges: $a = (3/32)(3\,\text{Tr}\,R^3 - \text{Tr}\,R)$, $c = (1/32)(9\,\text{Tr}\,R^3 - 5\,\text{Tr}\,R)$
- Hofman-Maldacena: $1/2 \leq a/c \leq 3/2$
