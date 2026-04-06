# Task 1.1: Simple Gauge Groups by Rank

## Classification of Simple Lie Algebras

Every simple Lie algebra over $\mathbb{C}$ belongs to one of the following families:

### Classical Series

| Cartan type | Group | Rank | Dimension | Dual Coxeter $h^\vee$ | Condition |
|-------------|-------|------|-----------|----------------------|-----------|
| $A_r$ | $\mathrm{SU}(r+1)$ | $r \geq 1$ | $r(r+2)$ | $r+1$ | |
| $B_r$ | $\mathrm{SO}(2r+1)$ | $r \geq 2$ | $r(2r+1)$ | $2r-1$ | $B_1 \cong A_1$ |
| $C_r$ | $\mathrm{Sp}(2r)$ | $r \geq 2$ | $r(2r+1)$ | $r+1$ | $C_1 \cong A_1$, $C_2 \cong B_2$ |
| $D_r$ | $\mathrm{SO}(2r)$ | $r \geq 3$ | $r(2r-1)$ | $2r-2$ | $D_3 \cong A_3$ |

### Exceptional Groups

| Cartan type | Rank | Dimension | Dual Coxeter $h^\vee$ |
|-------------|------|-----------|----------------------|
| $G_2$ | 2 | 14 | 4 |
| $F_4$ | 4 | 52 | 9 |
| $E_6$ | 6 | 78 | 12 |
| $E_7$ | 7 | 133 | 18 |
| $E_8$ | 8 | 248 | 30 |

### Low-Rank Isomorphisms

| Rank | Isomorphisms |
|------|-------------|
| 1 | $A_1 \cong B_1 \cong C_1$: all are $\mathrm{SU}(2) \cong \mathrm{SO}(3) \cong \mathrm{Sp}(2)$ |
| 2 | $B_2 \cong C_2$: $\mathrm{SO}(5) \cong \mathrm{Sp}(4)$ |
| 3 | $D_3 \cong A_3$: $\mathrm{SO}(6) \cong \mathrm{SU}(4)$ |

When these coincide as Lie algebras, they appear as a single entry. However, the global form of the group and representation content can differ (e.g. spinor reps of $\mathrm{SO}(5)$ vs fundamental of $\mathrm{Sp}(4)$), so both may be relevant for enumerating matter content.

## Groups by Rank

### Rank 1
- $A_1 = \mathrm{SU}(2)$, dim = 3, $h^\vee = 2$

### Rank 2
- $A_2 = \mathrm{SU}(3)$, dim = 8, $h^\vee = 3$
- $B_2 = \mathrm{SO}(5)$, dim = 10, $h^\vee = 3$
- $C_2 = \mathrm{Sp}(4)$, dim = 10, $h^\vee = 3$ (isomorphic to $B_2$ as Lie algebra)
- $G_2$, dim = 14, $h^\vee = 4$

Inequivalent algebras at rank 2: $A_2$, $B_2 \cong C_2$, $G_2$ (3 algebras, but $B_2$ and $C_2$ kept separate for representation enumeration).

### Rank 3
- $A_3 = \mathrm{SU}(4)$, dim = 15, $h^\vee = 4$
- $B_3 = \mathrm{SO}(7)$, dim = 21, $h^\vee = 5$
- $C_3 = \mathrm{Sp}(6)$, dim = 21, $h^\vee = 4$
- $D_3 = \mathrm{SO}(6)$, dim = 15, $h^\vee = 4$ (isomorphic to $A_3$)

Inequivalent: $A_3 \cong D_3$, $B_3$, $C_3$.

### Rank 4
- $A_4 = \mathrm{SU}(5)$, dim = 24, $h^\vee = 5$
- $B_4 = \mathrm{SO}(9)$, dim = 36, $h^\vee = 7$
- $C_4 = \mathrm{Sp}(8)$, dim = 36, $h^\vee = 5$
- $D_4 = \mathrm{SO}(8)$, dim = 28, $h^\vee = 6$
- $F_4$, dim = 52, $h^\vee = 9$

$D_4$ has triality symmetry (outer automorphism $S_3$): the vector **8**$_v$, spinor **8**$_s$, and conjugate spinor **8**$_c$ are permuted.

### Rank 5
- $A_5 = \mathrm{SU}(6)$, dim = 35, $h^\vee = 6$
- $B_5 = \mathrm{SO}(11)$, dim = 55, $h^\vee = 9$
- $C_5 = \mathrm{Sp}(10)$, dim = 55, $h^\vee = 6$
- $D_5 = \mathrm{SO}(10)$, dim = 45, $h^\vee = 8$

### Rank 6
- $A_6$, $B_6$, $C_6$, $D_6$
- $E_6$, dim = 78, $h^\vee = 12$

### Rank 7
- $A_7$, $B_7$, $C_7$, $D_7$
- $E_7$, dim = 133, $h^\vee = 18$

### Rank 8
- $A_8$, $B_8$, $C_8$, $D_8$
- $E_8$, dim = 248, $h^\vee = 30$

### Rank $\geq 9$
Classical series only: $A_r$, $B_r$, $C_r$, $D_r$.

## LiE Naming Convention

LiE uses Cartan labels directly:

| LiE name | Group |
|----------|-------|
| `A1` | $\mathrm{SU}(2)$ |
| `A2` | $\mathrm{SU}(3)$ |
| `B2` | $\mathrm{SO}(5)$ |
| `C2` | $\mathrm{Sp}(4)$ |
| `G2` | $G_2$ |
| `D4` | $\mathrm{SO}(8)$ |
| `F4` | $F_4$ |
| `E6` | $E_6$ |

This is the input format for `hilbert.py`'s `group()` function (currently hard-coded to `C`).

## Practical Scope

For computational feasibility, the project will proceed rank by rank. At each rank, all inequivalent simple algebras are processed. The number of algebras per rank:

| Rank | Classical | Exceptional | Total (inequivalent) |
|------|-----------|-------------|---------------------|
| 1 | 1 ($A_1$) | 0 | 1 |
| 2 | 2 ($A_2$, $B_2$) | 1 ($G_2$) | 3 |
| 3 | 3 ($A_3$, $B_3$, $C_3$) | 0 | 3 |
| 4 | 4 ($A_4$, $B_4$, $C_4$, $D_4$) | 1 ($F_4$) | 5 |
| 5 | 4 ($A_5$, $B_5$, $C_5$, $D_5$) | 0 | 4 |
| 6 | 4 | 1 ($E_6$) | 5 |
| 7 | 4 | 1 ($E_7$) | 5 |
| 8 | 4 | 1 ($E_8$) | 5 |
| $r \geq 9$ | 4 | 0 | 4 |

Note: at rank 2, $B_2 \cong C_2$ as Lie algebras but we list both for representation enumeration purposes. The "inequivalent" count above treats them as one algebra. In practice we enumerate representations for each Cartan type separately.
