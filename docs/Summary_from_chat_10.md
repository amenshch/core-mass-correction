# Handoff summary — RT mass-correction paper & pipeline

This conversation became very long. This document is a self-contained
summary to continue smoothly in a new chat. The project's persistent
memory (userMemories) already has earlier background; this summary
focuses on what changed in *this* session, since much of it is recent
and substantial.

## Files to upload in the new chat

At minimum, upload these so Claude can pick up exactly where this left off:

1. **`rt_mass_correction_method.tex`** — the paper, fully rewritten this
   session (4-axis method, new notation, 13 pages).
2. **`rt_mass_correction.bib`** — bibliography (includes Menshchikov2021,
   the getsf methodology paper, now cited for the "footprint factor" term).
3. **`mass_correction_pipeline_4d.py`** — the production pipeline,
   renamed and updated this session to compute the 4-axis correction
   alongside the original 3-axis one (see below). This is the file to
   keep developing, not the old `mass_correction_pipeline.py` name.
4. **`make_bes_grid.py`** — grid-construction script; docstring updated
   with a notation-mapping note (functionally unchanged, verified by diff).
5. **`make_flowchart_mpl.py`** — flowchart generator, fully updated to
   the 4-axis method and current notation.
6. Model grid catalogs (from the GitHub repo, `cats/` directory):
   - `bes_model_grid_final2_catalog` + `bes_model_grid_final2_recovery_tables_v2`
     (main production grid, 432 nodes / 8070 samples post pathology-fix)
   - `bes_model_grid_subdivision2_catalog` + `_recovery_tables`
     (the NEW diverse-concentration subdivision grid, 226 detectable
     models across 5 background columns — replaces the old 36-model
     critical-only strip)
7. **`loo4d_full_results.csv`** — the completed 432-node LOO comparison
   (3-axis vs. 4-axis), if further analysis of it is needed.
8. Figures if visual continuity matters: `fig_validation.pdf`,
   `fig_degeneracy.pdf`, `fig_factor_decomposition.pdf`,
   `fig_grid_coverage.pdf`, `fig_method_flowchart.pdf`.

The repo is `amenshch/Aquila-cutout-test` on GitHub if raw-fetching is
easier than uploading.

## What happened this session, in order

### 1. Found a real degeneracy in the 3-observable corrector
Near-neighbor test on the training grid: samples from *different* physical
models with near-identical (Σ_cloud, concentration, R_foot) can still need
true correction factors differing by 2.5x+ (e.g., two models with masses
0.105 and 2.109 M☉, same three observables, true C_M = 1.03 vs 2.67). The
spread does **not** shrink as neighbors get arbitrarily close — a genuine
degeneracy, not sampling/interpolation error. Traced to missing absolute
physical scale (the three observables are deliberately scale-invariant;
C_M is not, because it depends on self-shielding column ∝ ρ_c·r0).

### 2. Added a 4th observable: physical size
- Grid side: R_BE (physical truncation radius, parsecs) — already tabulated
  in the catalog (column 11, "R_BE_pc").
- Real-source side: **footprint radius** = √(FOOA·FOOB)/2, converted to
  parsecs at the source's **known distance**.
- Tested four candidate real-source definitions (half-max deconvolved,
  half-max beam-convolved, footprint radius, physical mass not tried
  yet) — footprint radius ties with deconvolved half-max as best, and is
  the more conceptually consistent match to R_BE (both are edge/extent
  measures, not half-max measures). **Adopted: footprint radius.**
- This reintroduces a genuine distance dependence, unlike the pure
  3-axis design — discussed explicitly in the paper (Appendix B extended
  to reconcile this with the earlier rejection of mass-indexed matching).

### 3. Found and mostly fixed a SED-fit robustness bug
Root cause: when per-band recovered flux fractions are highly imbalanced
across the 3 fitted bands (or when only 2 bands remain), the least-squares
modified-blackbody fit can be driven to an unphysical temperature at the
edge of the fitted grid, giving masses wrong by up to 800x. Fix: require
≥3 bands with recovered fraction ≥ MIN_FRAC (0.25, already an existing
constant reused for a new purpose) before a training sample is kept;
otherwise discard it. Reduced catastrophic (>10x) errors from 13/432 nodes
to 1/432. Cost: training corpus 8686→8070 samples, 432→406 nodes retain
≥1 usable sample. One residual case (node 211010) has moderate individual
fractions and isn't caught by any reasonable floor — documented as a known
limitation, not chased further.

### 4. New, better subdivision validation grid
Old one: 36 models, ALL at the critical truncation only (ξ_max=6.451) at
one background column — an unnoticed bug, only a 1-D slice of concentration
space was ever tested. New one (`bes_model_grid_subdivision2_catalog`):
226 detectable models spanning sub-critical through deep-super-critical,
at 5 background columns (geometric means between the 6 main-grid columns).
Generated via `make_bes_grid.py`'s own `build_node()`, unmodified — just
looped over new Σ_emb values — so it's exactly consistent with the main
grid's physics.

### 5. Full validation, 3-axis vs 4-axis, three independent tests
All numbers below use one representative sample per node where relevant
(closest R_foot to 1.84, the grid's median) to avoid pseudo-replication
from counting multiple truncation radii of the same physical model as
independent test cases.

| test | n | quantity | uncorrected | 3-axis | 4-axis |
|---|---|---|---|---|---|
| Injection (real getsf, real noise) | 108 | median | 0.76 | 0.97 | 0.96 |
| | | scatter (factor) | 1.81 | 1.93 | 1.81 |
| | | within 2× | 86% | 83% | 86% |
| Subdivision (new, diverse grid) | 166 | median | 0.84 | 1.02 | 1.01 |
| | | scatter (factor) | 1.56 | 1.58 | 1.32 |
| | | within 2× | 87% | 90% | 98% |
| Node-LOO (production grid) | 388 | median | 0.84 | 1.02 | 1.01 |
| | | scatter (factor) | 1.51 | 1.62 | 1.42 |
| | | within 2× | 87% | 88% | 95% |

Key finding stated plainly in the paper: the 3-axis corrector's scatter is
**worse than not correcting at all** in every test (median bias is fixed,
but per-source spread widens). The 4-axis corrector doesn't have this
problem — its scatter is at or below the uncorrected baseline in every
test. This was an important, non-obvious result worth remembering: adding
the 4th axis isn't just an incremental improvement, it changes the
3-axis method from "questionably beneficial" to "clearly beneficial" by
fixing something the median-only statistics were hiding.

### 6. Full notation overhaul (paper + flowchart)
| old | new | meaning |
|---|---|---|
| `concentration` (as the specific axis) | ζ (zeta) | peak-to-mean concentration |
| `R_foot` | φ (phi) | "footprint factor" (Men'shchikov 2021 term) |
| φ(ξ_max) [projected-width function] | renamed to `\varphi` | frees up φ for the line below |
| FWHM, AFWHM | H | half-maximum size |
| FOOA, FOOB | A_F, B_F | footprint major/minor axis |
| "SD" (surface-density band) | Σ | its common name |
| prof_λ(r), prof_SD(r) | I_λ(r), Σ(r) | actual profile names |
| RADMC-3D | *radmc-3d* (lowercase italic) | |
| "RT" | "radiative transfer" (expanded) | |
| literal beam values (13.5″, 18.2″, 24.9″, 36.3″) | O_λ, with O_160=13.5″ etc. | angular resolution symbol |

Also added missing formal definitions for F_flux, F_bg, F_temp, M_ref
(none were actually defined before, just used) — F_flux and F_bg have
exact formulas from the code (`M_SED3bs/M_rec` and `FXSDbsl/FXSDbs`
respectively); F_temp is the residual, F_temp ≡ C_M/(F_flux·F_bg).

**Important care point**: `\varphi` (the *original* projected-width
function, now renamed) must not be confused with the *new* `\phi`
(footprint factor, was R_foot) — these are deliberately different
symbols now. Also: the *scripts* (`make_bes_grid.py` etc.) still use the
*old* internal variable name `phi(xi)` for what the paper now calls
`\varphi` — a notation-mapping note was added to the script docstring to
prevent this exact confusion, but the actual variable name was **not**
renamed in code (too risky for a live pipeline; see below).

### 7. Pipeline integration (this was the last major piece of work)
`mass_correction_pipeline.py` → renamed **`mass_correction_pipeline_4d.py`**.

Before integrating, found and fixed a real naming collision: the code
already had an *unrelated*, experimental "4D" concept (`mass4d_interp`,
using the source's own reported mass as a 4th predictive feature — this
is actually the mass-indexed approach the paper's Appendix B explicitly
rejects). Renamed that to `mass_idx_interp`/`midx` with a clear docstring
explaining what it is, before introducing the real 4-axis corrector under
a non-colliding name.

What's now implemented in the production file:
- `build_recovery_samples()` extracts `rbe_pc` (R_BE in parsecs) per
  training sample.
- New `cm_size_interp(S)` builds the actual adopted 4-axis corrector,
  `(Σ, ζ, φ, R_BE) → C_M`.
- `load_getsf()` now returns `foot_radius_as` (footprint radius in
  arcsec, derived algebraically from the already-computed footprint
  area — no new column parsing needed: `sqrt(area/π)`).
- `correct_cloud()` computes **both** the 3-axis and 4-axis correction
  per source (converting footprint radius to pc using its own `distance`
  argument), and writes **both** to the output catalog as new columns
  `C_M_4D`, `M_CORR_4D` — **alongside**, not replacing, the original
  `C_M`/`M_CORR`. This was a deliberate choice: didn't want to silently
  change default behavior for anything downstream without an explicit
  decision from the user first.
- Verified end-to-end: ran the full 108-source injection test through
  the integrated `correct_cloud()` and got median=0.961, scatter=1.81,
  within-2×=86% — **exact match** to the standalone validation numbers
  in the table above. Not just "runs without crashing" — functionally
  verified correct.

**Not yet done** (flagged, not forgotten):
- The CLI `test` subcommand (`run_test()`, used for subdivision-grid
  validation) still only reports the 3-axis result. Natural next
  extension.
- Whether `C_M_4D`/`M_CORR_4D` should eventually become the *primary*
  `C_M`/`M_CORR` — not decided, left as an explicit open choice.
- Physical mass (as opposed to physical size) as an alternative/
  additional 4th-axis candidate was discussed as worth testing but never
  actually tried.
- The one residual degenerate case (node 211010, moderate-not-extreme
  band fractions) not caught by the MIN_FRAC floor.

## Known-good scripts for reproducing results (may need re-deriving if
context is lost, but were written and validated this session)
- Full 432-node LOO comparison script (checkpointed, ~50 min runtime,
  produces `loo4d_full_results.csv`) — was handed to the user directly
  as `loo_4axis_full.py`; if not saved locally, would need rewriting
  following the same method described in section 5 above.
- The near-neighbor degeneracy test (pairs samples from different nodes,
  bins by proximity, checks whether true C_M spread shrinks) — logic is
  described fully in the paper's Sect. "A fourth observable" and
  `fig_degeneracy.pdf`'s generation; would need rewriting if not
  otherwise available.

## Paper structure (current, 13 pages, compiles cleanly)
1. Introduction
2. The model grid (BE formalism, geometry justification, radiative
   transfer, grid coverage — largely unchanged from before this session
   except notation)
3. The mass correction
   - observables + recovery degeneracy (extended to flag the residual
     3-axis degeneracy)
   - the per-band footprint model (F_flux/F_bg/F_temp now properly defined)
   - interpolation and output constraints (SED-fit robustness fix
     documented here)
   - **NEW: "A fourth observable: physical size"** — the core new section
4. Validation (fully rewritten: 3 tests, 3-axis vs 4-axis table, the
   "why 3-axis doesn't reduce scatter but 4-axis does" discussion, the
   factor-decomposition and output-cap subsections kept but explicitly
   flagged as historical/superseded diagnostics that motivated the 4th
   axis rather than current recommendations)
5. Discussion (resolution dependence, low-mass edge — mostly unchanged)
6. Conclusions (updated for 4-axis)
Appendices: mass conservation, why-not-indexed-on-mass (extended),
interpolation accuracy, resolution-dependence estimate.

## Things worth double-checking at the start of the next chat
- Word count / page budget still within the 12+6 target (was 13 pages
  total at last compile, comfortably within budget).
- Whether the flowchart's validation numbers need yet another refresh if
  anything changes.
- Whether `run_test()` extension or the primary-output decision should
  be tackled next, or whether the user wants to return to the paper
  (e.g., a final full read-through) or start the companion application
  paper (seven-cloud results, deferred throughout).
