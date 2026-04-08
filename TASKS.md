# Task Tracker

## Section 0: Warm-Up — Review and Conventions

| Task | Description | Status | Output |
|------|-------------|--------|--------|
| 0.1 | Review a-maximization (hep-th/0304128) | Done | [task_0_1](tasks/task_0_1_a_maximization_review.md) |
| 0.2 | Review CMNS landscape paper (2408.02953) | Done | [task_0_2](tasks/task_0_2_CMNS_landscape_review.md) |
| 0.3 | Set notation and conventions | Done | [task_0_3](tasks/task_0_3_notation_conventions.md) |
| 0.4 | Document Nf=2N.nb Mathematica code | Done | [task_0_4](tasks/task_0_4_Nf2N_notebook.md) |
| 0.5 | Document hilbert.py | Done | [task_0_5](tasks/task_0_5_hilbert_py.md) |

## Section 1: Group Theory Data

| Task | Description | Status | Output |
|------|-------------|--------|--------|
| 1.1 | Enumerate simple gauge groups by rank | Done | [task_1_1](tasks/task_1_1_simple_groups.md) |
| 1.2 | Tabulate representation data (dim, T(R), reality) | Done | [task_1_2](tasks/task_1_2_rep_data.md), `data/all_reps.json` |
| 1.3 | GIO construction method (Hilbert series vs hand-coded) | Done | [task_1_3](tasks/task_1_3_gio_from_hilbert.md) |

## Section 2: Seed Theory Enumeration

| Task | Description | Status | Output |
|------|-------------|--------|--------|
| 2.1 | Asymptotic freedom: enumerate all matter content with b_0 > 0 | Done | `data/enumerate_seeds.py`, `data/seed_matters.json` |
| 2.2 | Anomaly cancellation: cubic (A_r) + Witten (Sp(N)) | Done | Integrated into `data/enumerate_seeds.py` |
| 2.3 | Filter for interacting IR fixed points via a-maximization | Done | `data/filter_seeds.py`, `data/consistent_seeds.json` |
| 2.4 | Catalog seed theories with operator spectra | Partial | `src/deform/analyze_seed.py` (needs orbit-aware PE) |

**Seed counts (2641 total):** A1:9, A2:17, A3:137, A4:81, B2:139, B3:199, C2:68, C3:117, D4:1841, G2:20, F4:13

## Section 3: Fixed-Point Analysis Tools

| Task | Description | Status | Output |
|------|-------------|--------|--------|
| 3.1 | Generalized a-maximization (FindCharges) | Done | `src/amax/FindCharges.wl` |
| 3.2 | Hilbert series for any simple group | Done | `src/hilbert/hilbert.py` |
| 3.3 | Identify relevant/super-relevant deformations | Partial | PE-based classification works; explicit monomial generation needs orbit-aware PE |
| 3.4 | Handle operator decoupling | Done | R <= 2/3 flagged, iteration stops |
| 3.5 | Limitations documented | Done | In [task_0_2](tasks/task_0_2_CMNS_landscape_review.md) and README.md |

## Section 4: Deformation Procedure

| Task | Description | Status | Output |
|------|-------------|--------|--------|
| 4.1 | Direct deformation (dW = O) | Done | `src/deform/iterate.py` |
| 4.2 | Flip deformation (dW = M*O) | Done | `src/deform/iterate.py` |
| 4.3 | Consistency checks | Done | In `src/amax/FindCharges.wl` |
| 4.4 | Equivalence via graph isomorphism | Done | `src/deform/equivalence.py` |

## Section 5: Iterative Tree Construction

| Task | Description | Status | Output |
|------|-------------|--------|--------|
| 5.1 | Depth-by-depth iteration engine | Done | `src/deform/iterate.py` |
| 5.2 | Parallel a-maximization | Done | multiprocessing.Pool, all CPUs |
| 5.3 | Track RG flow graph | Not started | |

**Test runs:**
- A1 (SU(2)): 7 non-decoupled seeds iterated to depth 5. Total 243 theories.
- G2: 20 seeds → 10 consistent at depth 1 (with dedup).

## Section 6: Database and Output

| Task | Description | Status | Output |
|------|-------------|--------|--------|
| 6.1 | SQLite database | Not started | |
| 6.2 | Validation against CMNS | Not started | |
| 6.3 | Summary tables and visualization | Not started | |

## Critical Fix Needed

### Orbit-Aware PE for Explicit Monomial Generation

**Problem:** Current pipeline has two broken approaches:
1. PE degree-based (`iterate.py`): gives counts, not explicit monomials. Misses singlets, M*meson, quartic.
2. Hand-coded (`OperatorSpectrum.wl`): SU(2)-specific. Misses baryons for SU(N>=3).

**Solution:** Orbit-aware PE (planned, not yet implemented)
1. Group fields into orbits under residual flavor symmetry
2. Run per-type PE with orbits as types → (order+1)^N_orbits cost
3. Use LiE `alt_tensor`/`sym_tensor` to find invariant tensor powers (group-independent)
4. Expand PE degrees into explicit monomials via combinatorial selection
5. F-term reduce via GroebnerBasis

**New files needed:**
- `src/hilbert/invariants.py` — find primitive invariants via LiE
- `src/hilbert/monomials.py` — expand PE degrees to explicit monomials
- Modified `src/deform/iterate.py` — orbit decomposition
- Modified `src/amax/OperatorSpectrum.wl` — accept external monomial list

**Validation:** SU(2) 8f at depth 1: must produce 134 relevant operators (27 mesons + M + M^2 + 15 M*meson + 90 quartic), matching OperatorSpectrum.wl output.
