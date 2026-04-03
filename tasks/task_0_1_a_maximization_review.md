# Task 0.1: Review of $a$-Maximization (hep-th/0304128)

**Paper:** K. Intriligator and B. Wecht, "The Exact Superconformal R-Symmetry Maximizes $a$"

## Core Problem

Every 4d N=1 SCFT has a unique superconformal U(1)_R inside SU(2,2|1). When the theory has flavor symmetries F_I, the R-symmetry is not uniquely fixed by anomaly cancellation alone — any linear combination

$$R_t = R_0 + \sum_I s_I F_I$$

is a valid candidate. The paper provides a universal principle to determine the correct superconformal R-symmetry.

## The $a$-Maximization Principle

**Theorem:** The superconformal U(1)_R is the one that **locally maximizes** the trial $a$-function:

$$a_{\text{trial}}(s) = \frac{3}{32}\left(3\,\text{Tr}\,R_t^3 - \text{Tr}\,R_t\right)$$

This is a **cubic polynomial** in the mixing parameters $s_I$.

### Stationarity Conditions

$$\frac{\partial a}{\partial s_I} = 0 \implies 9\,\text{Tr}(R^2 F_I) = \text{Tr}(F_I) \quad \forall I$$

### Local Maximum Condition

$$\frac{\partial^2 a}{\partial s_I \partial s_K} = \frac{27}{16}\,\text{Tr}(R\,F_I\,F_K) < 0 \quad \text{(negative definite)}$$

The negative definiteness follows from unitarity: the 2-point function coefficient $\tau_{IK} \propto -\text{Tr}(R\,F_I\,F_K)$ must be positive definite.

## Central Charges from R-Charges

For a gauge theory with gauge group $G$ (dimension $|G|$) and chiral multiplets $\Phi_i$ in representations $R_i$ (dimension $|R_i|$) with R-charges $r_i$:

$$\text{Tr}\,R^3 = |G| + \sum_i N_i(r_i - 1)^3 |R_i|$$

$$\text{Tr}\,R = |G| + \sum_i N_i(r_i - 1) |R_i|$$

The gaugino has R = 1 (contributes $(1-1)^3 = 0$ to Tr R^3 but is counted in |G|). Matter fermions have R-charge $(r_i - 1)$ where $r_i$ is the scalar component's R-charge.

**Central charges:**

$$a = \frac{3}{32}(3\,\text{Tr}\,R^3 - \text{Tr}\,R)$$

$$c = \frac{1}{32}(9\,\text{Tr}\,R^3 - 5\,\text{Tr}\,R)$$

## ABJ Anomaly Cancellation

For the R-symmetry to be non-anomalous under gauge interactions:

$$T(\text{adj}) + \sum_i N_i(r_i - 1)\,T(R_i) = 0$$

where $T(R)$ is the Dynkin index. This is equivalent to the vanishing of the NSVZ exact beta function at the fixed point, and **eliminates one R-charge variable**.

## Superpotential Constraints

Each superpotential term $W_k$ must have R-charge exactly 2:

$$R[W_k] = 2$$

Each such constraint is **linear** in the $r_i$ and eliminates one more variable. The superpotential also breaks some flavor symmetries, reducing the dimension of the mixing parameter space.

## Unitarity Bound

For gauge-invariant chiral primary operators:

$$R[\mathcal{O}] \geq \frac{2}{3} \quad \Longleftrightarrow \quad \Delta[\mathcal{O}] \geq 1$$

Saturation ($R = 2/3$, $\Delta = 1$) means the operator is a **free field**.

**Important:** This constrains gauge-invariant composites, NOT elementary field R-charges directly.

### When the Bound is Violated

If $a$-maximization yields $R[\mathcal{O}] < 2/3$ for some GIO, the operator must decouple as a free field. An accidental U(1) symmetry emerges, and one must redo $a$-maximization with the constraint $R[\mathcal{O}] = 2/3$.

## Connection to the $a$-Theorem

The $a$-maximization principle provides a near-proof of $a_{\text{IR}} < a_{\text{UV}}$:
- UV: $a_{\text{trial}}$ is maximized over the full flavor symmetry space.
- A relevant deformation breaks some flavor symmetries, restricting the parameter space.
- IR: $a_{\text{trial}}$ is maximized over a **smaller** subspace.
- Maximum over a subspace $\leq$ maximum over the full space.

Caveats: accidental symmetries in the IR can enlarge the space; only a local maximum is guaranteed.

## Algebraicity

Since $a$-maximization reduces to solving polynomial equations with rational coefficients, all superconformal R-charges are **algebraic numbers**. Consequently, $a$, $c$, and all scaling dimensions $\Delta = (3/2)R$ are algebraic.

## Free Field Check

For free chiral multiplets, extremizing $a_{\text{trial}}$ gives $(r_i - 1)^2 = 1/9$, i.e. $r_i = 2/3$ (local max) or $r_i = 4/3$ (local min). The correct free-field R-charge is $r = 2/3$, confirming maximization.

$$a_{\text{free}} = \frac{3}{16}|G| + \frac{1}{48}M$$

where $M$ is the number of free chiral multiplets.

## Key Formulas for Our Project

| Quantity | Formula |
|----------|---------|
| Trial $a$ | $\frac{3}{32}(3\,\text{Tr}\,R^3 - \text{Tr}\,R)$ |
| Central charge $c$ | $\frac{1}{32}(9\,\text{Tr}\,R^3 - 5\,\text{Tr}\,R)$ |
| ABJ constraint | $T(\text{adj}) + \sum_i N_i(r_i-1)T(R_i) = 0$ |
| Superpotential | $R[W_k] = 2$ for each term |
| Unitarity bound | $R[\mathcal{O}] \geq 2/3$ |
| Dimension relation | $\Delta = \frac{3}{2}R$ (chiral primaries) |
| Relevant operator | $R[\mathcal{O}] < 2$ |
| Free field | $R = 2/3$, $\Delta = 1$ |
