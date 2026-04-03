# Task 0.2: Review of Large Landscape of 4d SCFTs (2408.02953)

**Paper:** M. Cho, K. Maruyoshi, E. Nardoni, J. Song, "Large Landscape of 4d Superconformal Field Theories from Small Gauge Theories"

## Overview

Systematic classification of 4d N=1 SCFTs starting from rank-1 and rank-2 Lagrangian gauge theories. By iteratively applying relevant deformations and flippings, the authors identify **7,346 inequivalent fixed points** cataloged at `https://qft.kaist.ac.kr/landscape`.

## Seed Theories

Starting points: Lagrangian gauge theories at conformal fixed points with W = 0.

| Seed | Gauge Group | Matter |
|------|-------------|--------|
| 1 | SU(2) | 1 adjoint (**3**) + 2 fundamentals (**2**) |
| 2 | SU(2) | 1 adjoint (**3**) + 4 fundamentals (**2**) |
| 3 | SU(3) | 1 adjoint (**8**) + 1 fund (**3**) + 1 anti-fund (**3-bar**) |
| 4 | SO(5) | 1 antisymmetric (**10**) + 1 vector (**5**) |
| 5 | SO(5) | 1 antisymmetric (**10**) + 2 vectors (**5**) |
| 6 | SO(5) | 1 traceless symmetric (**14**) + 1 vector (**5**) |
| 7 | SO(5) | 1 traceless symmetric (**14**) + 2 vectors (**5**) |
| 8 | Sp(2) | 1 irrep [0,2] (**14**) + 2 fundamentals (**4**) |
| 9 | Sp(2) | 1 adjoint (**10**) + 2 fundamentals (**4**) |
| 10 | G_2 | 1 adjoint (**14**) + 1 fundamental (**7**) |

These are chosen to be inside the conformal window. All have simple gauge groups of rank 1 or 2.
**Every seed includes fundamental/vector representations** alongside a higher-dimensional representation (adjoint, antisymmetric, or symmetric tensor).
Note: SO(5) ~ Sp(2) as Lie algebras, but listed separately due to distinct representation content (**5** vs **4**, **10** vs **10** in different branches).

## The Deformation Algorithm

### Step 1: Operator Enumeration

At each fixed point, enumerate all gauge-invariant chiral operators (GIOs) and compute their R-charges using a-maximization.

Classify:
- **Relevant:** R[O] < 2 (can appear as dW = O)
- **Super-relevant:** R[O] < 4/3 (can be flipped: dW = M * O)

### Step 2a: Direct Deformation

Given a fixed point with superpotential W, for each relevant GIO O with R[O] < 2:
- Add dW = O to the superpotential: W -> W + O
- The new superpotential term imposes a new constraint: R[O] = 2
  (since every term in W must have R-charge 2)
- This constraint is linear in the trial R-charges and eliminates one free parameter
  from the a-maximization
- Redo a-maximization with this additional constraint to obtain new R-charges and new (a, c)
- The new constraint also breaks some flavor symmetries (those under which O is charged),
  reducing the dimension of the trial R-symmetry space

**Effect on the chiral ring:** the new superpotential term generates new F-term relations
d(W + O)/dPhi_i = 0 for each elementary field Phi_i.
These may set some operators to zero in the chiral ring that were nonzero before.

### Step 2b: Flip Deformation

**What flipping is:** Given a GIO O, flipping means introducing a new gauge-singlet
chiral field M and adding the superpotential coupling dW = M * O.
M is a **free field in the UV** with R\_UV[M] = 2/3 (the free-field R-charge).
The coupling M * O has R\_UV = 2/3 + R[O] and is relevant (R < 2)
precisely when R[O] < 4/3. At the IR fixed point, M acquires an anomalous dimension
and has R\_IR[M] = 2 - R\_IR[O].

**What flipping does to the chiral ring:**
- M is an elementary field, so its F-term equation is dW/dM = O = 0.
  This sets the composite operator O to **zero** in the chiral ring.
  O is removed from the spectrum entirely — it is not an operator of the new theory.
- The elementary fields composing O still exist, but their specific combination
  forming O is now a chiral ring relation.
- M itself is a new operator in the chiral ring with R[M] = 2 - R[O].
- F-terms of other elementary fields Phi_i gain new contributions M * dO/dPhi_i,
  which may generate additional chiral ring relations.

**Central charges:** In general, redo a-maximization with the new field M and new constraint.
When O is neutral under all residual flavor symmetries (so the flip does not change other R-charges),
the shift is purely from adding the singlet M:

> da = (3/32)(1 - R[O])(3R[O]^2 - 6R[O] + 2)
>
> dc = (1/32)(1 - R[O])(9R[O]^2 - 18R[O] + 4)

**Derivation:** M contributes d(Tr R^3) = (R[M]-1)^3 = (1-R[O])^3 and
d(Tr R) = (R[M]-1) = (1-R[O]) to the anomaly traces. Substituting into
da = (3/32)(3 d(Tr R^3) - d(Tr R)) gives the formula.
This is exact only when no other R-charges change; otherwise, redo full a-maximization.

### Step 3: Operator Decoupling

After a-maximization, if any GIO O\_d has R[O\_d] < 2/3,
the unitarity bound is violated and O\_d must decouple.

**Physical picture:**
- An accidental U(1) symmetry emerges in the IR under which O\_d is charged.
- The true superconformal R-symmetry mixes with this accidental U(1) to correct R[O\_d] to exactly 2/3.
- O\_d becomes a free field (Delta = 1) decoupled from the interacting sector.
- The remaining interacting sector has its own (corrected) R-charges and central charges.

**Operational implementation (CMNS):**
Decoupling is implemented by **flipping** O\_d: add a singlet X (free field at UV, R\_UV[X] = 2/3) with dW = X * O\_d.
- The coupling X * O\_d has R\_UV = 2/3 + R[O\_d] < 4/3 < 2, so it is relevant.
- F-term dW/dX = 0 sets O\_d = 0 in the chiral ring, removing it from the interacting sector.
- At the IR: R\_IR[X] = 2 - R[O\_d], with 4/3 < R[X] < 2 (relevant but not super-relevant).
- Redo a-maximization with the new superpotential.
- In CMNS, this does **not** count as increasing depth (it is a consistency correction, not a chosen deformation).

**Our project:** if any operator decouples (R < 2/3), we mark the theory as `"operator decoupled"` and **stop** — we do not flip or continue iterating from that theory. This matches the `Nf=2N.nb` implementation.

### Step 4: Iterate

Apply Steps 1-3 to every new fixed point. Stop when no new non-trivial fixed points are generated. The iteration depth defines "generations."

## Consistency Checks (Section 2.2)

Every candidate fixed point is tested against:

1. **Unitarity:** R[O] >= 2/3 for all chiral primaries (after decoupling)
2. **Hofman-Maldacena:** 1/2 <= a/c <= 3/2
3. **Positive central charges:** a, c > 0
4. **a-theorem:** a\_child < a\_parent

### Superconformal Index Tests

The paper uses the **superconformal index** for detailed checks:

> I(t,y;x) = Tr (-1)^F t^{3(R+2j1)} y^{2j2} x^f

**Reduced index** (strips derivative descendants):

> I\_hat = (1 - t^3 y)(1 - t^3/y)(I - 1)

From the reduced index:
- **Relevant operators:** coefficient of t^{3R} (with R < 2) at y^0
- **Conformal manifold dimension:** at order t^6, the coefficient counts (marginal operators) - (conserved currents). That is, alpha = n\_marginal - n\_currents. If alpha > 0, there is a conformal manifold.
- **Conserved currents:** coefficient of t^6 from adjoint of flavor symmetry

### N=2 Enhancement Detection

- **Sufficient:** positive coefficient of t^7(y + 1/y)
- **N=3:** coefficient 2 for t^7(y + 1/y)

**Critical limitation for our project:** We use the Hilbert series instead of the superconformal index. This has major consequences:

1. **Unitarity check is incomplete.** The superconformal index captures all short multiplets (including those built from derivatives and anti-chiral fields), while the Hilbert series only sees the chiral ring. An operator can satisfy the chiral ring unitarity bound R >= 2/3 but still signal pathology at the index level (e.g., higher-spin conserved currents indicating a free theory, or negative index coefficients). **This is a significant difference** — we may accept theories that the full index would reject, or miss decoupling that the index would detect.

2. **N=2/N=3 enhancement:** not detectable (requires t^7(y+1/y) coefficient).

3. **Conserved currents and conformal manifold dimension:** not computable (requires the full t^6 sector of the index, not just chiral ring).

4. **Emergent symmetries:** not visible — the index detects accidental symmetries via conserved current multiplets, which the Hilbert series cannot see.

## Equivalence Criteria

Two fixed points are **inequivalent** if they differ in:
- Central charges (a, c)
- Operator spectrum (sorted R-charges of GIOs)
- Superconformal index (when computed)

Theories reached by different paths can be identified if all invariants match.

**Non-commutativity:** Different orderings of the same set of deformations can lead to **inconsistent theories** (failing consistency checks at intermediate steps), even if some orderings succeed. That is, the path matters — a deformation that is consistent at one fixed point may not be consistent at another, so the landscape graph is sensitive to the order in which deformations are applied.

## Key Results

| Quantity | Value |
|----------|-------|
| Total fixed points | 7,346 |
| a/c range | 0.7228 to 1.2100 |
| Minimal SCFT | (a,c) = (633/2000, 683/2000) |
| Theories with a = c | Found (holographic relevance) |
| Conformal manifolds (alpha > 0) | Found |
| Possible new N=2 SCFT | 1 candidate |

## What Our Project Generalizes

| Aspect | CMNS | Our Project |
|--------|------|-------------|
| Gauge groups | 5 specific (SU(2), SU(3), SO(5), Sp(2), G_2) | All simple groups |
| Rank | 1-2 | Arbitrary |
| Matter | Small, specific reps | All reps with T(R) <= T\_max |
| GIO enumeration | Superconformal index | Hilbert series + F-terms |
| Unitarity check | Full (via index: all short multiplets) | **Partial** (chiral ring only — major limitation) |
| Emergent symmetry | Detectable via index | Not detectable |
| N=2/N=3 enhancement | Detectable | Not detectable |
| Conformal manifold | Computable (alpha from index) | Not computable |
| Seed theories | 10 hand-picked | Systematically enumerated |

## Notation Adopted from CMNS

- R-charges: R[O], with Delta = (3/2)R
- Relevant: R < 2; super-relevant: R < 4/3
- Flip: dW = M * O, with R[M] = 2 - R[O]
- Depth = number of superpotential terms
- Central charges: a = (3/32)(3 Tr R^3 - Tr R), c = (1/32)(9 Tr R^3 - 5 Tr R)
- Hofman-Maldacena: 1/2 <= a/c <= 3/2
