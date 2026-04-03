# Landscape: Classification of 4d N=1 SCFTs

## Project Overview

Systematic classification of 4d N=1 supersymmetric gauge theories with simple gauge group that flow to non-trivial interacting IR fixed points. This generalizes [Cho-Maruyoshi-Nardoni-Song (2408.02953)](https://arxiv.org/abs/2408.02953) to all simple gauge groups with all allowed matter representations.

**Method:** Start from seed theories (W=0), enumerate relevant deformations via Hilbert series + F-term reduction, deform by a single superpotential term, iterate up to depth 5.

## Repository Structure

```
Landscape/
├── CLAUDE.md              # This file
├── tasks.tex / tasks.pdf  # Research plan and task list
├── references/            # Reference papers and code
│   ├── Nf=2N.nb           # Mathematica: a-maximization + deformation pipeline (SU(N))
│   ├── hilbert.py         # Python: Hilbert series via LiE + plethystic logarithm
│   ├── 2408.02953.pdf     # CMNS landscape paper
│   └── hep-th0304128.pdf  # Intriligator-Wecht a-maximization paper
├── src/                   # Source code (generalized pipeline)
│   ├── hilbert/           # Hilbert series computation (LiE + Mathematica)
│   ├── amax/              # a-maximization solver
│   ├── deform/            # Deformation engine (iterate depth 0→5)
│   └── db/                # Database interface (SQLite)
├── data/                  # Group theory data and representation tables
├── results/               # Output databases and exports per gauge group
└── notebooks/             # Analysis notebooks (Mathematica / Jupyter)
```

## Key Conventions (following 2408.02953)

- **Representations:** specified by Dynkin labels `[λ_1, ..., λ_r]`
- **R-charges:** from a-maximization; Δ = (3/2)R for chiral primaries
- **Unitarity bound:** R[O] >= 2/3
- **Relevant deformation:** R[O] < 2
- **Super-relevant / flippable:** R[O] < 4/3
- **Hofman-Maldacena bound:** 1/2 <= a/c <= 3/2
- **Depth:** number of superpotential terms

## Pipeline Overview

1. **Hilbert series** (LiE + Python multiprocessing): PE computes gauge-invariant content for any simple group + Dynkin-label representations. Truncated at chosen order in fugacities.
2. **a-maximization** (Mathematica): ABJ constraint + W constraints → maximize a_trial. Exact Solve, fallback to NSolve with WorkingPrecision→50.
3. **F-term reduction** (Mathematica): GroebnerBasis of ∂W/∂Φ_i = 0, PolynomialReduce to reduce chiral ring.
4. **Operator classification:** Read all GIOs from PE with R < 2 (relevant) and R < 4/3 (flippable). Deduplicate via graph isomorphism (toGlobalGraph).
5. **Iteration:** depth 0 → depth 5, storing results in SQLite at each level.

## External Dependencies

- **LiE** (Lie group computations): used by hilbert.py for plethystic exponential
- **Mathematica / wolframscript**: used for a-maximization, plethystic logarithm, GroebnerBasis
- **Python 3**: multiprocessing for parallel PE computation

## Development Notes

- When modifying `hilbert.py`, the key generalization points are: `group()` (accept all Lie types, not just C_r) and `matter()` (accept arbitrary Dynkin labels, not 9 fixed slots).
- The `FindCharges` function in `Nf=2N.nb` is the template for a-maximization; generalize by parameterizing Dynkin indices μ_i and dimensions d_i from LiE output.
- Emergent symmetries and N=2/N=3 enhancement are NOT detectable without the full superconformal index. This is a known limitation.
- Theories are identified up to equivalence via superpotential graph isomorphism (toGlobalGraph).
