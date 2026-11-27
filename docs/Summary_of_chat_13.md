# Summary of chat 13 — for transition to chat 14

Written 2026-08-14. Everything below is measured, not assumed, unless flagged
otherwise. Numbers marked PENDING will change when the extended grid is
computed.

---

## 1. WHAT THIS CHAT ACHIEVED

The chat began with an unsolved problem: injected model cores did not resemble
real Aquila cores, which prevented any meaningful injection test. That problem
is solved. Along the way, five defects in the pipeline and the test machinery
were found and fixed, and the model grid was extended.

**Headline results**

| Quantity | At the start | Now |
|---|---|---|
| Overlap of injected and real populations in (column, contrast) | 0.318 | 0.890 |
| Contrast, Kolmogorov–Smirnov statistic vs real | — | 0.052 (p = 0.29) |
| Four-observable corrector, node leave-one-out scatter | 1.42 | 1.11 |
| Four-observable corrector, within a factor of 2 | 95% | 100% |
| Production grid | 720 nodes (518 built) | 4446 possible, 2443 built |

---

## 2. THE PLACEMENT METHOD (the chat's main product)

Five requirements, each established by measurement:

1. **Selection must not use the quantity under test.** 
   Matching injected cores
   to real ones by reported mass builds the bias under test into the test
   sample. Match instead on what the extraction measures well — background
   column PEAK^BGF, contrast PEAK^SBF/PEAK^BGF, and observed size AFWHM — and
   let mass follow. Mass enters the node score with weight 0.1 against 1.0 for
   the others; raising it to parity degrades the contrast match from KS 0.063
   to 0.174.

2. **The embedding window must tile the column axis.** 
   Grid rungs are spaced by
   a factor 2, so the window must be one octave wide and centred:
   `SD_emb/√2 ≤ Σ_local < √2·SD_emb`. The previous one-sided window was half a
   rung wide, leaving half of every octave unreachable; 520 of 1030 real
   sources, including the population median, fell in those gaps, and models
   drawn for them were silently reassigned upward. Fixing this raised
   reachability from 50% to 99%.

3. **Contrast must be scored as C − 1, not C.** C spans only 0.09 dex over the
   real sources; C − 1 = PEAK^SRC/PEAK^BGF spans 0.87 dex.

4. **Sizes and separations must use observables, not R_BE.** 
   R_BE reaches 428″
   for models whose observed FWHM is 25–50″. Selecting or separating on R_BE
   discarded every model above 1.8 M☉ and gave survivors enormous exclusion
   zones. Verified that model FWHMSDbs and getsf AFWHM are the same quantity:
   over 451 matched cores the ratio has median 1.001, scatter 0.127 dex, valid
   below ~50″ (above which the recovered size saturates near 50–60″ through
   footprint truncation).

5. **Each core must be steered to its own drawn column.** Steering all cores of
   one node to their median collapses the within-node spread and halves the
   contrast–column correlation.

**Separation rule.** Exclusion radius = the radius enclosing 95% of the model's
own flux (from `bes_model_grid_final2_profiles`), never less than the observed
footprint 0.92·FWHMSDbs. Neither geometric radius works: node i03j18k08 has
R_BE = 631″ but 95% of its flux within 43″. Applies to every pair, and the same
radius (scaled) governs clearance from real source footprints.

**Clustering.** `MAX_FOOTPRINT_OVERLAP = 0.30`, chosen to match the real
sources: 44% of real Aquila sub-field sources overlap their nearest neighbour's
footprint, 22% by more than 30%. Forbidding overlap entirely makes the injected
sample cleaner than the data (5th-percentile separation 1.91 vs real 0.78 in
units of the pair's mean FWHM).

---

## 3. THE FOUR-WAY INJECTION DESIGN

One placement rendered four times. Every core appears at the same position with
the same model in all four; they differ only in the two switches.

| Run | Cloud fluctuations | Neighbouring cores | Config |
|---|---|---|---|
| as placed | present | present | `inject_config_v3.py` |
| flattened | removed | present | `inject_config_v3_flat.py` |
| isolated | present | removed | `inject_config_v3_iso.py` |
| flattened+isolated | removed | removed | `inject_config_v3_flat_iso.py` |

`CONTROL_ISOLATED` filters the completed placement to a maximal
non-overlapping subset, so the isolated truth table is a strict subset of its
companion's. Runs are independent — separate directories are correct; pairing
comes from the shared seed, not shared files.

**Verification before analysing any set:**
```
diff <(grep -v '^#' inj_..._iso_truth.txt | awk '{print $2,$3,$4}' | sort) \
     <(grep -v '^#' inj_..._truth.txt     | awk '{print $2,$3,$4}' | sort)
```
Only `>` lines are acceptable.

**Critical lesson:** compare core by core, never sample against sample. An
earlier flattened pair compared sample-wise appeared to show a 39% bias from
cloud fluctuations; it vanished when the same cores were compared in both.

---

## 4. VALIDATION RESULTS (pipeline v6, SNR 1.5, MIN_SED_BANDS 3)

All entries restricted to sources inside the correction domain. Ratio =
mass obtained / reference mass; scatter = std of log10(ratio) as a factor.

| Test | n | Treatment | Median | Scatter | Within 2× |
|---|---|---|---|---|---|
| Node leave-one-out | 388 | uncorrected | 0.838 | 1.36 | 89% |
| | | 3-observable | 1.130 | 1.35 | 95% |
| | | **4-observable** | **1.013** | **1.11** | **100%** |
| Subdivision grid | 156 | uncorrected | 0.838 | 1.33 | 90% |
| | | 3-observable | 1.168 | 1.34 | 96% |
| | | **4-observable** | **1.024** | **1.11** | **100%** |
| Background scan | 110 | uncorrected | 0.756 | 1.98 | 79% |
| | | 3-observable | 1.087 | 1.94 | 88% |
| | | 4-observable | 0.934 | 1.95 | 85% |
| Population injection (flat+iso, 4 seeds) | 95 | uncorrected | 0.79 | 1.56 | 85% |
| | | 4-observable | 0.928 | 1.56 | 88% |

The cSD injection test is **retired** — superseded grid, superseded placement,
blocks 04/05 break the extraction.

### Contamination, measured core by core (4 seeds pooled)

| Comparison | Pairs | Scatter 1 | Scatter 2 | Shift (dex) | p |
|---|---|---|---|---|---|
| Fluctuations, neighbours present | 113 | 1.84 | 1.58 | −0.006 | 0.99 |
| Fluctuations, neighbours removed | 76 | 1.96 | 1.46 | +0.029 | 0.07 |
| Blending, fluctuations present | 108 | 2.34 | 2.40 | −0.001 | 0.02 |
| Blending, fluctuations removed | 95 | 1.57 | 1.56 | −0.009 | 0.0008 |

Cloud fluctuations: scatter only (factor ~1.2–1.3), no bias.
Blending: ~2% mass suppression, established; interacts with fluctuations.

### Seven clouds (PENDING regeneration with extended grid)

| Cloud | d (pc) | sources | in domain | med C_M 3-obs | med C_M 4-obs | med M_ex (M☉) |
|---|---|---|---|---|---|---|
| Scorpius | 130 | 187 | 24 | 1.443 | 1.238 | 0.082 |
| Ophiuchus | 139 | 965 | 232 | 1.422 | 1.228 | 0.033 |
| Aquila | 260 | 990 | 853 | 1.436 | 1.242 | 0.131 |
| Orion A | 432 | 1445 | 775 | 1.425 | 1.227 | 0.675 |
| California | 470 | 991 | 257 | 1.421 | 1.262 | 0.605 |
| Cygnus X | 1400 | 1874 | 1558 | 1.408 | 1.193 | 3.399 |
| W3/W4/W5 | 1700 | 1654 | 587 | 1.427 | 1.256 | 4.722 |

C_M spans only 1.193–1.262 over 130–1700 pc; correlation with log d = −0.13.
**Distance invariance confirmed on real catalogues.**

Values are lower than the ~2.9–3.6 quoted before this chat, because two-band
fits (which attracted huge spurious corrections) are now excluded.

---

## 5. PIPELINE CHANGES (v5 → v6)

1. **`IN_DOMAIN` flag.** `_HullTest` (convex-hull facet inequalities on
   standardised axes) reports whether a query interpolates or extrapolates.
   This restores a filter that was lost silently when the estimator changed
   from Delaunay (returned NaN outside the hull) to LocalLinear (extrapolates).
   Leave-one-out admits exactly 388 of 403 nodes — the same count Delaunay gave.
   New output columns: `FOOT_AS, IN_HULL_3D, IN_HULL_4D, IN_DOMAIN,
   SUPPORT_3D, SUPPORT_4D`.
   NB: `LL_BOUND`/`EXTRAP_*` is a *different* test (did the plane leave its
   neighbours' value range) and fires almost never; it is not a domain test.

2. **`MIN_SED_BANDS = 3`** (was hard-coded 2). A modified blackbody has two
   free parameters, so a two-band fit has zero degrees of freedom, passes
   exactly through both points, and chi-squared is zero by construction —
   `QUALITY = ok` cannot catch it. **Every** core wrong by >5× was fitted on
   two bands; none fitted on three or four was.

3. **`SNR_MIN_PEAK = SNR_MIN_TOTAL = 1.5`** (was 2.0), chosen by sweep:

| S/N | bands | cores | median | scatter | wrong >5× |
|---|---|---|---|---|---|
| 1.0 | 2 | 427 | 1.003 | 3.87 | 13% |
| 1.0 | 3 | 331 | 0.906 | 2.30 | 6% |
| 1.5 | **3** | **230** | **0.815** | **1.77** | **1%** |
| 2.0 | 3 | 135 | 0.751 | 1.52 | 0% |

Note the trap: at the loosest setting the median looks perfect (1.003) but is
large errors of both signs cancelling.

---

## 6. THE EXTENDED GRID

`bes_grid_final3.txt`: 2443 models, 2087 detectable, 396 hole.
Ladders: SD_emb 9 values 3.75e20–9.6e22 (×2); M 26 values 0.0016–30.6 M☉
(×1.483); FWHM 19 values 0.003–0.25 pc (×1.279) = 2.4–198″ at 260 pc.
Contrast floor 1.05 (lowered from 1.1: only 1.0% of real sources lie below
1.05 vs 6.2% below 1.10).

**Gap closure, faintest contrast available vs real requirement (~1.23–1.26):**

| SD_emb | 15–26″ | 26–35″ | 35–45″ | 45–55″ | 55–80″ | usable |
|---|---|---|---|---|---|---|
| 3.75e20 | 3.92 | 4.10 | 3.02 | 2.99 | 2.13 | no |
| 7.5e20 | 2.46 | 2.55 | 2.01 | 2.00 | 1.56 | no |
| 1.5e21 | 1.73 | 1.77 | 1.50 | 1.50 | 1.28 | marginal |
| 3.0e21 | 1.31 | 1.38 | 1.25 | 1.25 | 1.14 | yes |
| 6.0e21 | 1.16 | 1.13 | 1.13 | 1.12 | 1.07 | yes (was 1.87!) |
| ≥1.2e22 | ~1.06 | ~1.07 | ~1.06 | ~1.06 | ~1.05 | yes |

**HARD PHYSICAL LIMIT FOUND.** The mass ladder cannot go lower: at every size
the minimum mass sits at ξ_max ≈ 0.50–0.55, the bottom of the tabulation, which
is the uniform-sphere limit. A BE sphere of given FWHM has a minimum mass
(M ∝ c_s²R/G). Consequently **a Bonnor–Ebert sphere cannot look like a real
low-contrast core on a very faint background** — at 3.75e20 the faintest
possible 26–35″ sphere still shows contrast 4.1 where real sources show 1.25.
This is a limitation of the model family, to be stated in the paper, not a grid
defect. Worth one check: whether tabulating ξ_max below 0.5 admits more nodes.

**Extrapolation diagnosis.** Of 8106 sources in seven clouds, 3820 (47%) fall
outside the training hull; 3287 fail on background column alone, ≤99 on any
other axis. The three new low rungs address this, but see the limit above.

**Hole nodes.** 13% of corrector training nodes, 15% of training samples — they
ARE used in every cloud correction (I initially claimed otherwise; wrong).
Tested: leave-one-out gives 1.012/1.11 with them, 1.001/1.12 without. Harmless,
keep them. Use `bsl` stamps, never `bs`.

---

## 7. SUBDIVISION GRID — MODIFICATIONS NEEDED

Must be recomputed: it is the geometric midpoints of the production ladders,
and all three production ladders changed.

In `make_bes_grid.py`:
```python
SD_VALUES = [5.30e20, 1.06e21, 2.12e21, 4.24e21,
             8.49e21, 1.70e22, 3.39e22, 6.79e22]       # 8 midpoints
M_MIN, M_MAX, N_M    = 0.00196, 25.12, 25             # 25 midpoints
FW_MIN, FW_MAX, N_FW = 0.003376, 0.2211, 18           # 18 midpoints
```
Product 8 × 25 × 18 = 3600; expect ~2000 to survive the ξ_max inversion and
the contrast floor (production: 2443 of 4446).
Compute this **last**, after the production ladders are final.

---

## 8. SCRIPT INVENTORY (proposed subdirectory layout)

### `method/` — the mass-correction method itself
| Script | Purpose |
|---|---|
| `mass_correction_pipeline_6.py` | **The pipeline.** Correctors, SED fit, hull test, `correct_cloud` |
| `fit_mass.py` | `CMInvariantCorrector` class |
| `getsf_columns.py` | Header-position column reader (avoids column-number fragility) |
| `getsf_bg.py` | Background helper used by `inject.py` |
| `read_getsf_sed.py` | SED reader |

### `grid/` — model grid creation
| Script | Purpose |
|---|---|
| `make_bes_grid.py` | **Grid definition.** Writes `bes_grid_final3.txt` |
| `grid_extension.py` | Computes/justifies the extended ladders, reports node counts |
| `flag_hole_models.py` | Recomputes the ρ_edge < ρ_emb flag from a catalogue |
| `tabulate_recovery.py` | Builds recovery tables |
| `add_recoverable_mass_multiband.py` | Recoverable-mass columns |
| `calibrate_floors.py` | Floor calibration |

### `inject/` — injection experiments
| Script | Purpose |
|---|---|
| `run_inject_v3.py` | **The injection driver.** Do not edit; all settings in configs |
| `inject.py` | Placement primitives (`valid_placement_mask`, `local_sigma`) |
| `inject_config_v3.py` | As placed |
| `inject_config_v3_flat.py` | Flattened background |
| `inject_config_v3_iso.py` | Isolated control |
| `inject_config_v3_flat_iso.py` | Both |
| `build_joint_target_v3.py` | Builds the real-source target table |
| `joint_target_Aquila_v3.txt` | 1030 real Aquila sources (data) |
| `dryrun_placement_v3.py` | Selection + placement, no images written |
| `smoke_test_v3.py` | End-to-end test of `main()` on synthetic inputs |

### `tests/` — validation and analysis
| Script | Purpose |
|---|---|
| `run_validation.py` | **All four validation tests**, scored twice (all / in-domain) |
| `analyse_quad_injection.py` | Four-way injection analysis, pools seeds |
| `paired_injection_v3.py` | Paired flattening comparison |
| `analyse_injection_v3.py` | Single-run injection analysis |
| `check_injection_match_v3.py` | Population match vs real sources |
| `sweep_sed_bands.py` | SED band-threshold sweep against known truth |
| `compare_size_measures.py` | Model FWHMSDbs vs getsf AFWHM calibration |
| `analyse_aquila_sim.py` | The original Könyves-era simulated field |
| `loo_4axis_full_v2.py` | Standalone leave-one-out (fast, runs locally) |
| `run_validation_tests.py` | Chat-10 driver — **superseded** by `run_validation.py` |

### `figures/` — paper figures
| Script | Produces |
|---|---|
| `aa_plotstyle.py` | Shared A&A conventions (`clean_log`, `no_top_minor`) |
| `make_fig_grid_coverage3.py` | Fig. 2, grid coverage (new grid) |
| `make_fig_grid_coverage.py` | previous version — superseded |
| `make_flowchart_mpl.py` | Fig. 1, method flowchart |
| `make_fig_footprint_model.py` | footprint model |
| `make_fig_degeneracy.py` | degeneracy |
| `make_fig_factor_decomposition.py` | factor decomposition |
| `make_fig_validation.py` | validation figure |
| `fig_sim_vs_real.py` | simulated vs real populations |
| `make_fig_capping_tradeoff.py` | capping trade-off |

### `misc/`
`cmf_mle.py`, `check_separation.py`, `estimator_compare.py`,
`make_compare_pseudoreplication.py`, `patch_*.py` (historical patches),
`reload_project.sh`

---

## 9. TODO LIST

**Immediate (next chat)**
1. Regenerate recovery tables for `bes_grid_final3.txt`.
2. Re-run all four validation tests; compare v6-old-grid vs v6-new-grid.
3. Re-run seven clouds; measure how much of the 47% extrapolation is recovered.
4. Compute the subdivision grid (Sect. 7 above) and re-run that test.
5. Update in the method paper: Table 1 (validation), Table 3 (contamination),
   Sect. 4.2 domain figures, Fig. 2. All carry `\draftnote{}` markers.
6. Update the seven-cloud table in the living document.

**Open questions**
- Does tabulating ξ_max below 0.5 admit further low-mass nodes?
- The three under-estimated outliers (ratio 0.04–0.16) are all large models
  43–70″ — footprint truncation, or the same photometry failure? Needs more
  seeds.
- Injected cores pass `QUALITY = ok` at 29% vs 46% for real sources. Why?
- Selection function: rank correlation of contrast on column is +0.51 in real
  sources but ~0 in recovered injected ones. Detection destroys it.
- Profile-fit fourth axis: 0.086 dex vs 0.045 dex now achieved by the
  footprint radius under the new estimator. The case for it has weakened;
  `profile_fit_test.py` was never preserved.

**Deferred**
- Forward-modelling the selection function to recover the true population
  (stated as beyond scope in the paper).
- Shortening the method paper to 12 + 6 pages (currently 18 total: main text
  through references ~14, appendices ~4).

---

## 10. CONVENTIONS AND GOTCHAS

**Papers.** `rt_mass_correction_method.tex` = the A&A method paper (has
"_method"). `rt_mass_correction.tex` = the living development document.

**LaTeX style** (from `svdensities.tex`): `\rm` not `\mathrm`; powers not
radicals, `(a+b)^{1/2}`; one-line formulas, two per numbered equation where
related; numbered sub(sub)sections, no run-in headings; italic only for
telescope/method names; panel labels as `(\emph{left})`; en-dashes `--`.
Macros: `\CM` = 𝒞_M, `\Fflux \Ftemp \Fbg` = ℱ_F ℱ_T ℱ_B, `\Mex \Mbe \Msed
\Mref \Mcorr \Mrec`, `H_{\rm BE}` (renamed from `H_intr`).
Keep `sec:`/`fig:`/`tab:`/`eq:`/`app:` label prefixes.

**Band indices differ between catalogues** — this caused two false alarms:
- Guoyin whole-field: 160, **161=SD**, 250, 255, 350, 500 → SD is band **02**
- Injection/sub-field: 070, 160, **161=SD**, 250, 350, 500 → SD is band **03**
- Aquila-sim: 070, 160, 165, 250, 350, 500 → SD is band **03**
Always read the `! wavelength` block before indexing.

**Stamps: use `bsl`, never `bs`.** `bs` subtracts a flat interpolated
background and goes negative inside R_BE for the 15% of nodes whose envelope is
rarer than the cloud. `bsl` subtracts the true crater floor and is non-negative.
The setting also selects the catalogue columns used for selection
(`ICSD<suffix>`, `M_SED4<suffix>`) so model chosen = model injected.

**`SET_TAG` is ignored** in core-mass-function mode; output names are built as
`cmf<slope>_n<cores>_s<seed>` with `_flat`/`_iso` appended automatically.

**GitHub directory prefixes are dates** and are not guessable
(260811, 260812, 260813, 260814…). List the directory via the web page when the
API rate limit is exhausted:
`curl -s https://github.com/amenshch/Aquila-cutout-test/tree/main/<dir>`

**Stale `/tmp/fit_mass.py`** can shadow `python/fit_mass.py`; delete it if
`CMInvariantCorrector` import fails.

**Column counts are 1-based** in conversation.


RT code radmc3d takes on input rho_c, xi_max, T_BE, and SD_emb, and they 
seem to have enough digits, then the running computation is safe, and nothing 
needs recomputing.

The four inputs are printed as:

Input	Format	Example	Significant digits
rho_c	{:>11.4e}	5.1892e-18	5
xi_max	{:>8.4f}	0.5508	4
T_BE	{:>7.3f}	11.608	5
SD_emb	"%.2e"	3.75e+20	3, but exact
SD_emb looks the coarsest but is not: the values are exactly the ladder 
entries 3.75×10²⁰, 7.5×10²⁰, 1.5×10²¹ and so on, all of which are representable 
to three digits, so nothing is lost.

The weakest is xi_max at four digits, which is a relative precision of 
about 2×10⁻⁴ at the low end. That is far below the RT noise and irrelevant 
to the resulting density structure.

So every model now being computed is the model it should be. The only defect 
is that the catalogue's M_BE and FWHM_pc columns are rounded labels on correctly 
computed spheres.

What this changes

It becomes purely a validation-accuracy issue rather than a physics one, 
and the fix is a reprint after the computation finishes:

Let the run complete.
Change ROW_FMT: M_BE from {:>10.4f} to {:>12.6g}, FWHM_pc from {:>8.4f}
to {:>10.6g}, and widen the matching HEADER_FMT fields.
Re-run make_bes_grid.py alone. It is analytic and takes seconds.
Verify the node count and every (i, j, k) are unchanged.

The reprinted M_BE is then the exact ladder value rather than a value 
rounded to ±3%, which matters because M_BE is the denominator of every 
validation ratio and the corrector is being quoted at the 1% level.

One consequence worth noting for the next chat: the existing 
bes_model_grid_final2_catalog has the same rounded columns, so a small part 
of the leave-one-out scatter quoted throughout this work — 1.11 for the 
four-observable corrector — is rounding in the truth rather than error in 
the method. It will be interesting to see whether that figure improves 
slightly with the reprinted grid, and it is a clean, independent reason to 
expect the new grid to validate better.
