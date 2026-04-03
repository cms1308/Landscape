# Task 0.3: Notation and Conventions

Following the conventions of 2408.02953 (CMNS).

## Gauge Group

- $G$: simple Lie group, specified by Lie type + rank (e.g. $A_4 = \mathrm{SU}(5)$, $C_2 = \mathrm{Sp}(2)$, $G_2$)
- $h^\vee$: dual Coxeter number
- $|G|$: dimension of $G$ (= dimension of adjoint)
- $T(R)$: Dynkin index of representation $R$, normalized so $T(\mathrm{fund}) = 1/2$ for $\mathrm{SU}(N)$
- $|R|$: dimension of representation $R$
- $A(R)$: cubic anomaly coefficient

## Representations

Specified by **Dynkin labels** $[\lambda_1, \ldots, \lambda_r]$ where $r = \mathrm{rank}(G)$. This is the input format for LiE.

## Matter Content

- Chiral multiplets $\Phi_i^{(a)}$ ($a = 1, \ldots, N_i$) in representation $R_i$ with R-charge $r_i$
- Gauge-singlet (flip) fields: $X_j$, $M_j$
- Multiplicity vector: $n = (N_{X}, N_{M}, N_{R_1}, N_{R_2}, \ldots)$

## R-charges and Scaling Dimensions

- $R[\mathcal{O}]$: superconformal R-charge of operator $\mathcal{O}$
- $\Delta[\mathcal{O}] = \frac{3}{2} R[\mathcal{O}]$ for chiral primaries
- Gaugino: $R = 1$
- Free chiral field: $R = 2/3$, $\Delta = 1$

## Central Charges

$$a = \frac{3}{32}(3\,\mathrm{Tr}\,R^3 - \mathrm{Tr}\,R)$$

$$c = \frac{1}{32}(9\,\mathrm{Tr}\,R^3 - 5\,\mathrm{Tr}\,R)$$

where

$$\mathrm{Tr}\,R^3 = |G| + \sum_i N_i (r_i - 1)^3 |R_i|$$

$$\mathrm{Tr}\,R = |G| + \sum_i N_i (r_i - 1) |R_i|$$

## Constraints

- **ABJ anomaly cancellation:** $T(\mathrm{adj}) + \sum_i N_i (r_i - 1) T(R_i) = 0$
- **Superpotential:** $R[W_k] = 2$ for each term $W_k$
- **Asymptotic freedom:** $b_0 = 3\,T(\mathrm{adj}) - \sum_i N_i\,T(R_i) > 0$

## Operator Classification

| Type | Condition | Role |
|------|-----------|------|
| Relevant | $R[\mathcal{O}] < 2$ | Direct deformation: $\delta W = \mathcal{O}$ |
| Super-relevant | $R[\mathcal{O}] < 4/3$ | Flippable: $\delta W = M \cdot \mathcal{O}$ |
| Marginal | $R[\mathcal{O}] = 2$ | Exactly marginal if on conformal manifold |
| Irrelevant | $R[\mathcal{O}] > 2$ | Cannot deform |
| Unitarity-violating | $R[\mathcal{O}] < 2/3$ | Must decouple |

## Bounds

- **Unitarity:** $R[\mathcal{O}] \geq 2/3$ for all gauge-invariant chiral primaries
- **Hofman-Maldacena:** $1/2 \leq a/c \leq 3/2$
- **$a$-theorem:** $a_\mathrm{IR} < a_\mathrm{UV}$ along any RG flow

## Terminology

- **Depth** (= length): number of superpotential terms. Seed theories have depth 0.
- **Flip:** introduce free singlet $M$ with $\delta W = M \cdot \mathcal{O}$. Sets $\mathcal{O} = 0$ in chiral ring via F-term.
- **Seed theory:** asymptotically free gauge theory with $W = 0$ flowing to interacting IR fixed point.
- **GIO:** gauge-invariant operator (element of the chiral ring).
- **Chiral ring:** ring of gauge-invariant chiral operators modulo F-term relations $\partial W / \partial \Phi_i = 0$.
