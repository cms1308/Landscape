# Landscape

Systematic classification of 4d N=1 superconformal field theories with simple gauge group, generalizing [Cho-Maruyoshi-Nardoni-Song (2408.02953)](https://arxiv.org/abs/2408.02953) to all simple gauge groups with all allowed matter representations.

## Method

Starting from seed theories (asymptotically free gauge theories with no superpotential), we:

1. Perform **a-maximization** to determine superconformal R-charges at the IR fixed point
2. Enumerate **gauge-invariant operators** via the Hilbert series (plethystic exponential)
3. Classify operators as relevant (R < 2), super-relevant (R < 4/3), or unitarity-violating (R <= 2/3)
4. **Deform** by a single superpotential term (direct or flip) and repeat
5. Iterate up to **depth 5** (five superpotential terms)

## Current Status

### Completed

| Component | Description | Location |
|-----------|-------------|----------|
| Group enumeration | All simple Lie groups by rank | `data/enumerate_reps.py`, `data/all_reps.json` |
| Representation data | Dynkin labels, T(R), dim, reality for all reps with T < 3h^v | `data/all_reps.json` |
| Seed enumeration | All AF matter content with anomaly cancellation (cubic + Witten) | `data/enumerate_seeds.py`, `data/seed_matters.json` |
| Seed filtering | a-maximization on all candidates | `data/filter_seeds.py`, `data/consistent_seeds.json` |
| a-maximization | Generalized FindCharges for any simple group | `src/amax/FindCharges.wl` |
| Hilbert series | Generalized PE computation via LiE | `src/hilbert/hilbert.py` |
| Operator spectrum | GIO enumeration + F-term reduction | `src/amax/OperatorSpectrum.wl` |
| Deformation engine | Iterative depth-by-depth pipeline with parallel a-max | `src/deform/iterate.py` |
| Equivalence checking | Graph isomorphism for theory deduplication | `src/deform/equivalence.py` |

### Consistent Seeds (2641 total)

| Group | h^v | Candidates | Consistent |
|-------|-----|-----------|------------|
| A1 = SU(2) | 2 | 12 | 9 |
| A2 = SU(3) | 3 | 22 | 17 |
| A3 = SU(4) | 4 | 157 | 137 |
| A4 = SU(5) | 5 | 92 | 81 |
| B2 = SO(5) | 3 | 161 | 139 |
| B3 = SO(7) | 5 | 227 | 199 |
| C2 = Sp(4) | 3 | 80 | 68 |
| C3 = Sp(6) | 4 | 130 | 117 |
| D4 = SO(8) | 6 | 1961 | 1841 |
| G2 | 4 | 26 | 20 |
| F4 | 9 | 17 | 13 |

### In Progress

- Orbit-aware PE for explicit monomial generation (group-independent GIO construction)
- Full iteration to depth 5 for all groups

## Repository Structure

```
Landscape/
├── README.md
├── CLAUDE.md                 # Development context and conventions
├── tasks.tex / tasks.pdf     # Research plan
├── tasks/                    # Task writeups (reviews, documentation)
├── references/               # Reference papers and code
│   ├── Nf=2N.nb              # Reference: a-max + deformation pipeline (SU(N))
│   ├── hilbert.py            # Reference: Hilbert series via LiE
│   ├── 2408.02953.pdf        # CMNS landscape paper
│   └── hep-th0304128.pdf     # Intriligator-Wecht a-maximization paper
├── src/
│   ├── amax/                 # a-maximization (Mathematica)
│   │   ├── FindCharges.wl    # Generalized a-maximization solver
│   │   └── OperatorSpectrum.wl  # GIO enumeration + F-term reduction
│   ├── hilbert/              # Hilbert series (Python + LiE)
│   │   └── hilbert.py        # PE computation for any simple group
│   ├── deform/               # Deformation engine (Python + Mathematica)
│   │   ├── iterate.py        # Depth-by-depth iteration with parallel a-max
│   │   ├── equivalence.py    # Graph isomorphism deduplication
│   │   └── analyze_seed.py   # Seed analysis (a-max + PE + classification)
│   └── db/                   # Database (planned)
├── data/                     # Computed data
│   ├── all_reps.json         # Representation data for all groups
│   ├── seed_matters.json     # All AF matter content
│   ├── consistent_seeds.json # Seeds passing a-maximization
│   ├── enumerate_reps.py     # Script: enumerate reps with T < 3h^v
│   ├── enumerate_seeds.py    # Script: enumerate AF matter content
│   └── filter_seeds.py       # Script: filter seeds via a-max
├── results/                  # Output (per gauge group)
└── notebooks/                # Analysis notebooks
```

## Dependencies

- **Python 3** with `networkx`, `multiprocessing`
- **LiE** ([Lie group computations](http://www-math.univ-poitiers.fr/~maavl/LiE/))
- **Mathematica / wolframscript** (a-maximization, GroebnerBasis, plethystic logarithm)

## Usage

### Enumerate representations for a gauge group

```bash
python3 data/enumerate_reps.py G2 A2 B3
```

### Enumerate asymptotically free matter content

```bash
python3 data/enumerate_seeds.py G2 A2
```

### Filter seeds via a-maximization

```bash
python3 data/filter_seeds.py G2
```

### Run deformation iteration

```bash
# Group, max depth, [optional: single seed index]
python3 src/deform/iterate.py G2 5
python3 src/deform/iterate.py A1 5 1    # single seed
```

## Conventions

Following [Intriligator-Wecht (hep-th/0304128)](https://arxiv.org/abs/hep-th/0304128) and [CMNS (2408.02953)](https://arxiv.org/abs/2408.02953):

- **Representations** specified by Dynkin labels (LiE format)
- **R-charges** from a-maximization; scaling dimension Delta = (3/2)R
- **Unitarity bound:** R[O] >= 2/3 (saturation = free field, decoupled)
- **Relevant deformation:** R[O] < 2; add dW = O
- **Super-relevant / flippable:** R[O] < 4/3; add free singlet M with dW = M*O
- **Flip:** F-term dW/dM = O = 0 removes O from chiral ring; M is free at UV (R=2/3), acquires anomalous dimension at IR
- **Hofman-Maldacena bound:** 1/2 <= a/c <= 3/2
- **Depth:** number of superpotential terms
- **Decoupling:** if any GIO has R <= 2/3, theory is flagged and iteration stops (we do not flip decoupled operators)

## Limitations

- **No emergent symmetry detection:** the Hilbert series cannot detect accidental symmetries visible in the superconformal index
- **No N=2/N=3 enhancement detection:** requires the full superconformal index
- **Unitarity check is partial:** only checks chiral ring operators, not all short multiplets
- **Conformal manifold dimension not computable** without the index

## References

1. K. Intriligator and B. Wecht, "The exact superconformal R-symmetry maximizes a," [hep-th/0304128](https://arxiv.org/abs/hep-th/0304128)
2. M. Cho, K. Maruyoshi, E. Nardoni, and J. Song, "Large landscape of 4d superconformal field theories from small gauge theories," [2408.02953](https://arxiv.org/abs/2408.02953)
