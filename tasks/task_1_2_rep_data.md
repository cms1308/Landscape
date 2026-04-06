# Task 1.2: Representation Data

Representation data for simple Lie groups relevant to the landscape classification.
Dynkin index $T(R)$ normalized so that $T(\text{fund}) = 1/2$ for $\mathrm{SU}(N)$.

All data can be computed via LiE (`dim`, `contragr`, `alt_tensor`) but the Dynkin index normalization requires care — LiE's `norm()` uses a group-dependent inner product. The tables below use the standard physics convention.

## $A_1 = \mathrm{SU}(2)$, dim = 3, $h^\vee = 2$

| Dynkin label | Name | dim | $T(R)$ | Reality |
|-------------|------|-----|--------|---------|
| [1] | fund (**2**) | 2 | 1/2 | pseudo-real |
| [2] | adj (**3**) | 3 | 2 | real |
| [3] | **4** | 4 | 5 | pseudo-real |
| [4] | **5** | 5 | 10 | real |

## $A_2 = \mathrm{SU}(3)$, dim = 8, $h^\vee = 3$

| Dynkin label | Name | dim | $T(R)$ | Reality |
|-------------|------|-----|--------|---------|
| [1,0] | fund (**3**) | 3 | 1/2 | complex |
| [0,1] | anti-fund ($\bar{\mathbf{3}}$) | 3 | 1/2 | complex |
| [1,1] | adj (**8**) | 8 | 3 | real |
| [2,0] | sym (**6**) | 6 | 5/2 | complex |
| [0,2] | sym-bar ($\bar{\mathbf{6}}$) | 6 | 5/2 | complex |
| [3,0] | **10** | 10 | 15/2 | complex |

## $A_3 = \mathrm{SU}(4)$, dim = 15, $h^\vee = 4$

| Dynkin label | Name | dim | $T(R)$ | Reality |
|-------------|------|-----|--------|---------|
| [1,0,0] | fund (**4**) | 4 | 1/2 | complex |
| [0,0,1] | anti-fund ($\bar{\mathbf{4}}$) | 4 | 1/2 | complex |
| [0,1,0] | anti-sym (**6**) | 6 | 1 | real |
| [1,0,1] | adj (**15**) | 15 | 4 | real |
| [2,0,0] | sym (**10**) | 10 | 3 | complex |
| [0,0,2] | sym-bar ($\overline{\mathbf{10}}$) | 10 | 3 | complex |

## $B_2 = \mathrm{SO}(5)$, dim = 10, $h^\vee = 3$

| Dynkin label | Name | dim | $T(R)$ | Reality |
|-------------|------|-----|--------|---------|
| [1,0] | vector (**5**) | 5 | 1 | real |
| [0,1] | spinor (**4**) | 4 | 1/2 | pseudo-real |
| [0,2] | antisymmetric traceless (**10**) | 10 | 3 | real |
| [2,0] | symmetric traceless (**14**) | 14 | 5 | real |

Note: The adjoint of SO(5) is the antisymmetric tensor **10** = [0,2] with $T = h^\vee = 3$.

## $C_2 = \mathrm{Sp}(4)$, dim = 10, $h^\vee = 3$

| Dynkin label | Name | dim | $T(R)$ | Reality |
|-------------|------|-----|--------|---------|
| [1,0] | fund (**4**) | 4 | 1/2 | pseudo-real |
| [0,1] | traceless antisym (**5**) | 5 | 1 | real |
| [2,0] | adj = sym (**10**) | 10 | 3 | real |
| [0,2] | **14** | 14 | 5 | real |
| [1,1] | **16** | 16 | 7/2 | pseudo-real |

Note: $B_2 \cong C_2$. The vector **5** of SO(5) = [0,1] of Sp(4), and the spinor **4** of SO(5) = fund [1,0] of Sp(4).

## $G_2$, dim = 14, $h^\vee = 4$

| Dynkin label | Name | dim | $T(R)$ | Reality |
|-------------|------|-----|--------|---------|
| [1,0] | fund (**7**) | 7 | 1 | real |
| [0,1] | adj (**14**) | 14 | 4 | real |
| [2,0] | **27** | 27 | 9 | real |
| [1,1] | **64** | 64 | 28 | real |
| [0,2] | **77** | 77 | 42 | real |

Note: All representations of $G_2$ are real (the group has no complex or pseudo-real representations).

## $B_3 = \mathrm{SO}(7)$, dim = 21, $h^\vee = 5$

| Dynkin label | Name | dim | $T(R)$ | Reality |
|-------------|------|-----|--------|---------|
| [1,0,0] | vector (**7**) | 7 | 1 | real |
| [0,0,1] | spinor (**8**) | 8 | 1 | real |
| [0,1,0] | adj (**21**) | 21 | 5 | real |
| [2,0,0] | sym traceless (**27**) | 27 | 9 | real |

## $C_3 = \mathrm{Sp}(6)$, dim = 21, $h^\vee = 4$

| Dynkin label | Name | dim | $T(R)$ | Reality |
|-------------|------|-----|--------|---------|
| [1,0,0] | fund (**6**) | 6 | 1/2 | pseudo-real |
| [0,1,0] | traceless antisym (**14**) | 14 | 5/2 | real |
| [2,0,0] | adj = sym (**21**) | 21 | 4 | real |
| [0,0,1] | **14'** | 14 | 5/2 | pseudo-real |

## $D_4 = \mathrm{SO}(8)$, dim = 28, $h^\vee = 6$

| Dynkin label | Name | dim | $T(R)$ | Reality |
|-------------|------|-----|--------|---------|
| [1,0,0,0] | vector (**8**$_v$) | 8 | 1 | real |
| [0,0,1,0] | spinor (**8**$_s$) | 8 | 1 | real |
| [0,0,0,1] | conj-spinor (**8**$_c$) | 8 | 1 | real |
| [0,1,0,0] | adj (**28**) | 28 | 6 | real |

$D_4$ has triality: the three **8**-dimensional representations (vector, spinor, conj-spinor) are related by the $S_3$ outer automorphism and all have $T = 1$.

## $F_4$, dim = 52, $h^\vee = 9$

| Dynkin label | Name | dim | $T(R)$ | Reality |
|-------------|------|-----|--------|---------|
| [0,0,0,1] | fund (**26**) | 26 | 3 | real |
| [1,0,0,0] | adj (**52**) | 52 | 9 | real |
| [0,0,1,0] | **273** | 273 | 63 | real |

Note: All representations of $F_4$ are real.

## General Formulas

For $\mathrm{SU}(N)$:
- $T(\text{fund}) = 1/2$
- $T(\text{adj}) = N = h^\vee$
- $T(\text{sym}) = (N+2)/2$
- $T(\text{anti}) = (N-2)/2$

For $\mathrm{SO}(N)$:
- $T(\text{vector}) = 1$
- $T(\text{adj}) = N-2 = h^\vee$
- $T(\text{spinor})$: depends on $N$

For $\mathrm{Sp}(2r)$:
- $T(\text{fund}) = 1/2$
- $T(\text{adj}) = r+1 = h^\vee$

## How to Compute via LiE

```python
# Dimension
echo 'setdefault A3; print(dim([1,0,1]))' | lie

# Contragredient (conjugate representation)
echo 'setdefault A3; print(contragr([1,0,0]))' | lie

# Reality: check contragr(x)==x, then alt_tensor(2,x)|null(rank) for pseudo-real
echo 'setdefault A1; x=[1]; a=alt_tensor(2,x); print(a|null(Lie_rank))' | lie

# Adjoint highest weight
echo 'setdefault G2; print(adjoint())' | lie
```

The Dynkin index requires careful normalization (LiE's `norm()` uses a non-standard inner product for some Lie types). Use the formula $T(R) = h^\vee \cdot \text{dim}(R) \cdot C_2^{\text{LiE}}(R) / (\text{dim}(G) \cdot C_2^{\text{LiE}}(\text{adj}))$ with verification against known values.
