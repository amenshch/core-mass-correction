# Summary of chat 14

## 1. Validation of the mass correction, extended grid

Pipeline: `mass_correction_pipeline_6.py`. Corrector unchanged in form: local
linear regression on four observables — background surface density, concentration,
footprint factor, footprint radius in arcsec.

**Extended grid (final3), leave-one-node-out, 1597 nodes, 30 022 samples:**
four observables give median 1.010, scatter 1.13, all nodes within 2x, typical
error 0.027 dex; three observables 1.054 and 0.066 dex. Restricted to the volume
final2 spanned, 0.024 dex.

**Subdivision grid (1400 nodes, absent from training):** median 1.009, scatter
1.12, all within 2x, 0.028 dex. Interpolation between grid points is sound.

**Neighbour count:** over a factor of 32 (80 to 2560) the four-observable error
moves only 0.026 to 0.030 dex. The final2/final3 difference is not a smoothing
artefact: matching the neighbour *fraction* makes it slightly larger.

**Hole nodes:** four observables 0.037 dex against 0.173 for three; uncorrected
median 0.483. The angular truncation radius is what makes them correctable.

**Three-observable drift:** its median moves 0.04–0.08 between grids and 0.06
across the neighbour scan; the four-observable median moves under 0.012.

## 2. Two defects found and fixed in the pipeline

**Signal-to-noise threshold.** The code had drifted to 2.0 while the paper
adopts 1.5. Set to 1.5. This changed every real-data result. With it,
`sweep_sed_bands.py` reproduces Table 2 exactly (230 cores, median 0.815,
scatter 1.77, 2 cores wrong by >5x), and the background-scan row of Table 1
reproduces exactly (110 sources, 0.756, 1.98, 79%).

**Early binding of the grid file names.** `build_recovery_samples` declared its
two file names as default arguments, evaluated at import; assigning `mc.CAT67`
later had no effect and failed *silently*. Signature changed to
`(cat_path=None, rectab_path=None)` with the constants read in the body.

**Paper error corrected.** Table 1, population-injection row: the
three-observable column had duplicated the four-observable values. Correct
values are median 1.13, scatter 1.60, 84% within 2x. That row is the
flattened-and-isolated arm restricted to the 95 cores also recovered in the
flattened run.

## 3. Seven clouds, final3 corrector, threshold 1.5

**Correction (chat 24):** the medians below are those of the THREE-observable
factor (column C_M of the version-6 catalogs).  The adopted four-observable
factor (column C_M_4D there, C_M in version 7) has medians 1.28, 1.27, 1.25,
1.26, 1.30, 1.21, 1.30 for the same clouds (1.21-1.30).  The CMF slopes, the
within-cloud spread, and the bound fractions of Sect. 4 were computed with the
four-observable factor and stand.

Median correction factor: 1.36, 1.40, 1.43, 1.42, 1.41, 1.39, 1.45 for
Scorpius, Ophiuchus, Aquila, Orion A, California, Cygnus X, W3/W4/W5 — a spread
of 7 per cent over 130 to 1700 pc with no distance trend. On final2 Scorpius was
an outlier at 1.19, caused by the low-column grid edge.

The threshold moves the factor by at most 0.02 but raises the number of
corrected sources by 6 to 39 per cent per cloud.

**Surface-density band is not the same band number in every catalogue:** band 02
for five clouds, 03 for California, 01 for Cygnus X. The pipeline identifies it
from the data; any hard-coded band number would be wrong for two of seven.

## 4. Core mass function and boundedness

**CMF slopes do not steepen.** Changes are −0.01 to −0.08, i.e. very slight
*flattening*. Any claim of steepening is withdrawn. Reason: within a cloud the
correction spans only 1.14 to 1.65 between its 10th and 90th percentiles, and is
anti-correlated with mass (−0.36 to −0.56) and with background column (−0.41 to
−0.68).

**Bound fraction rises** from 45.0 to 50.2 per cent over 5574 resolved cores,
287 crossing the threshold: 11 per cent more candidate prestellar cores and 23
per cent more mass in them. Absolute levels use a radius derived from the
footprint factor and should be recomputed from catalogue sizes.

## 5. Inversion of the four observables

Node withheld, 30 022 samples: true mass 0.114 dex, xi_max 0.177, density
contrast 0.265, central density 0.336. The last two are **not determined** and
must not be quoted. Predicting a *ratio* is always better: C_M 0.028 against
0.114 for the mass; C_alpha 0.043 against 0.095; C_T 0.003 against 0.006.

**Stability classification.** Conventional criterion (deconvolved size,
alpha<=2) scores 73.4% against a 59.8% baseline, calling 44 per cent of
sub-critical cores bound. Recalibrating the threshold to ~0.8 gives 88.8%.
Grid inversion gives 89.5%, a branch classifier 89.2% (no gain). **At the
footprint factors real catalogues reach (1.8–2.1) all converge to 81–85%.**

The criterion has a **ceiling of 82.0%** even with the model's own mass, radius
and temperature, because the Bonnor-Ebert mass is double-valued in xi_max.

## 6. Temperature: a bias of mass, not of temperature

The reported `F_temp` is a residual `C_M/(F_flux F_bg)`, clamped below at unity;
13.5 per cent of samples are raised by the clamp. It is not a measurement.

Measured directly: true mass over SED mass with the correct background is
**1.080**; the corresponding temperature ratio is **1.001**. So the fitted
temperature is unbiased in the median while the mass from it is 8 per cent low.
Correcting the temperature cannot help. C_T fails on injections for exactly this
reason: predicted 0.997 (0.980–1.021) against required 0.923 (0.770–1.107),
rank correlation −0.22.

The warm bias seen on injected sources (8.3 per cent as placed, 2.4 per cent
flattened) comes from **background subtraction**, not from internal temperature
structure.

## 7. Sizes: the strongest transferable result

Target: intrinsic width at half maximum. Node withheld, grid: measured width
undeconvolved 0.907 and 0.124 dex; **Gaussian deconvolution 0.755 and 0.151 dex**
— worse than not deconvolving; **grid inversion 0.980 and 0.049 dex**.

**It transfers.** On getsf-extracted injected sources the grid inversion is best
in all four arms, median 1.028–1.067 and 0.013–0.033 dex, against 0.852–0.871
and 0.062–0.071 dex for Gaussian deconvolution.

The advantage grows as sources become marginally resolved: at H/O = 1.0–1.3 and
500 um, 0.021 dex against 0.301.

**Important qualification:** almost all the accuracy comes from the
surface-density morphology (concentration, footprint factor, footprint radius),
not from the low-resolution width. With only the background and the band's own
width the error is 0.114–0.154 dex, no better than reporting the measured width.
The method is therefore not a deconvolution of a marginally resolved
measurement; it recovers the intrinsic size from the shape of a resolved one and
requires a high-resolution surface-density map.

## 8. Two model-side gains that did NOT transfer

**Per-band width ratios.** On the grid they raise stability classification from
82–86 to ~96.5 per cent at realistic phi and reduce C_M from 0.028 to 0.015 dex
and xi_max from 0.177 to 0.015. On getsf-extracted sources they *degrade*
everything. Cause: a 1–4 per cent systematic offset between model-side and
extraction-side widths is enough, because the discriminating information lies in
a narrow range of the ratios.

**The error budget explains the limit.** On injected sources the four-observable
corrector reaches 0.122–0.127 dex against 0.027 in model space: real-data error
is **extraction-dominated by a factor of five**. No refinement of the corrector
can help while that holds.

## 9. Correcting the hires surface density

Bias measured on model stamps: 8 per cent at 3.75e20, 21 at 6e21, 52 at 4.8e22,
57 at 9.6e22 cm^-2.

**Mechanism, demonstrated numerically.** Two layers of equal column at 12 and
18 K, fitted with one modified blackbody over the four Herschel bands, return 75
per cent of the true column; 10 and 20 K give 59 per cent; a cold 8 K interior
with an 18 K skin gives 36–54 per cent. A sky whose bands are computed from a
surface-density map through one assumed temperature per pixel has no such spread
and needs no correction — which is why the simulated sky cannot test this.

**Calibration.** One row per pixel of model stamps; answer C_N =
N_true/N_hires. Leave-one-model-out on 34 961 pixels from 21 models: **the hires
surface density alone** gives 0.0014–0.0022 dex below contrast 2 and ~0.012
above 5. The intermediate hires fits add little except between contrast 2 and 5.

**Geometry.** Six models at 6e21 with cloud radii spanning 1.08e4 to 1.07e7 AU:
0.0022 dex with the column alone, and at matched *apparent* column the six agree
to 1–6 per cent. A more compact cloud is read lower AND needs more correction —
the two move together along one curve. Spherical geometry therefore suffices;
cylinders would add an unconstrained parameter.

**A gate on agreement between hires' own fits does not exist:** where the fits
agree most closely the required correction is still 1.083 and 99.6 per cent of
pixels need more than 2 per cent.

**Correction curve wobbles:** 502 of 1999 steps run downhill, largest reversal
0.17 in C_N, producing visible terracing along isophotes. Smoothing plus forced
monotonicity now applied, but the real remedy is more models, especially above
2e22 where the curve is barely constrained.

## 10. The injection test, and what it revealed about hires

First attempts were void: the reference image was itself a hires product, and
separately the two runs had been made with different versions of an FFT
convolution, giving a spurious signal 159 per cent of the injected column. Both
diagnosed and fixed.

**With the fixed code, the decisive test.** Growing apertures around 92 isolated
injected sources, difference of injected and original maps:

| aperture | 13.5" map | 36.3" map |
|---|---|---|
| 12" | 0.591 | 0.498 |
| 24" | 0.703 | 0.620 |
| 48" | 0.950 | 0.898 |
| **60"** | **1.018** | **0.991** |

**hires conserves the column but redistributes it outward over ~60 arcsec.**
There is no photometric loss. Every alarming number from smaller apertures — 0.42
inside footprints, 0.20 in the brightest bin, negative sums around 3e20 — was an
aperture artefact.

**Consequence for the correction:** near a compact source the reconstructed value
in a single pixel bears no simple relation to the true column at that position,
so a per-pixel correction cannot behave well there however good the calibration.
Pixel-by-pixel remains appropriate in the diffuse cloud.

## 11. Open items

- **G_0 dependence, untested and a priority.** 
  The grid is at G_0 = 10; Aquila is
  nearer 1. A stronger field widens the temperature spread without changing the
  apparent column much, so unlike geometry it may NOT be absorbed. Ten paired
  models at G_0 = 1 were chosen: i01j27k10, i02j15k12, i03j27k10, i04j16k13,
  i05j17k10, i07j19k10, i08j23k12, i09j19k09, plus i05j09k01 (unresolved) and
  i02j25k09 (contrast 36). The absolute fitted temperature is now carried in the
  pixel table as the observable that would absorb it. This affects the mass
  correction too, since the seven clouds are not in comparable fields.
- Incompleteness correction: the one application for which the grid is necessary
  rather than merely better. Forward-model framing agreed: measure detection
  probability against true mass and background column, then test which input mass
  functions are consistent with the observed catalogue, rather than reweighting.
- Bound fractions to be recomputed from catalogue sizes.
- The method paper needs reframing around the broader results, not C_M alone.

## 12. Handoff: analysing the hires method itself

Proposed next direction: rather than correcting hires output, examine the method.
Both inputs are in the project: `getsf.tex` and the `hires` script.

Known behaviour to start from, established above:

1. Column is conserved but **redistributed outward over ~60 arcsec**; peaks are
   fainter and faint sources can drop below detection.
2. The reconstruction **overshoots then undershoots around sharp edges**; the
   depression around bright peaks is a known problem and appears in models with a
   large density and temperature jump at the core edge.
3. hires needs **a background under the source**; on source-only images it fails.
4. The surface density it returns is **systematically low** because one modified
   blackbody is fitted per line of sight — a separate defect from 1–3, and the one
   the C_N calibration addresses.

Questions worth putting first: whether the redistribution scale is set by the
500 um beam and can be reduced; whether the median filtering of the temperature
map causes the leakage and the peak suppression; and whether the resolution
improvement and the temperature fitting can be reformulated so the bias in 4 does
not arise, rather than being corrected afterwards.

## Terminology

"Surface density" for the quantity; "column" only where it names the material
along a line of sight within a unit cross-section. Applied throughout the living
document and all 31 scripts.
