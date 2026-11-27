# Correction of Herschel-derived prestellar core masses

Radiative-transfer-calibrated correction of systematic biases in masses of prestellar cores extracted by *getsf* from *Herschel* images.
The correction factor C_M = M_BE / M_reported is predicted from four distance-invariant observables measured in the 13.5 arcsec surface density image: background surface density (PEAK^BGF), peak-to-mean concentration, footprint factor (footprint size over half-maximum size), and angular footprint radius.
It is calibrated on a grid of RADMC-3D models of Bonnor-Ebert spheres embedded in uniform envelopes and validated by leave-one-node-out tests, an independent subdivision grid, and injections into real Aquila images.

Earlier development history (chats 1-14 and the hires work) is in the archive repository `amenshch/Aquila-cutout-test`.

## Layout

| Directory | Contents |
|---|---|
| method/ | The correction pipeline and its helpers |
| grid/ | Grid definition, recovery-table construction, RADMC-3D setup |
| inject/ | Population-matched injection of model cores into real images |
| tests/ | Validation tests and their analysis |
| figures/ | Scripts producing the paper figures, shared plot style |
| misc/ | Auxiliary analysis (core mass function fits, estimator comparison) |
| hicorrect/ | Surface density correction of hires images (Fortran reference and Python version) |
| data/grid/ | Grid node lists, model catalogs, recovery tables, radial profiles |
| data/grid_G0/ | Grid ladder over the radiation field G_0 = 0.98-483 |
| data/clouds/ | getsf catalogs of seven clouds; Koenyves et al. (2015) Aquila tables |
| data/injections/ | Catalogs and truth tables of the injection runs and background scans |
| latex/ | Method paper, living document, bibliography, A&A style files, figures |
| docs/ | Chat summaries and working notes |

## Files

### method/
- `mass_correction_pipeline_6.py` — the pipeline: SED refit (S/N > 1.5 in peak and total flux, at least 3 bands), three- and four-observable correctors (local linear regression, K = 320), convex-hull domain flag `IN_DOMAIN`, per-cloud corrected catalog. Modes: `correct`, `loo`, `cmloo`, `test`.
- `fit_mass.py` — class `CMInvariantCorrector`.
- `getsf_columns.py` — reads getsf catalog columns by header position.
- `getsf_bg.py` — background helper used by the injection code.
- `read_getsf_sed.py` — reads getsf SED catalogs.

### grid/
- `make_bes_grid3.py` — defines the production grid; writes `bes_grid_final3.txt` (2443 models).
- `make_bes_grid3_subdivision.py` — defines the subdivision grid at the geometric midpoints of the production ladders.
- `tabulate_recovery.py` — builds the recovery tables (per node, per truncation radius, per band) from the model images.
- `flag_hole_models.py` — flags models whose boundary density is below the envelope density.
- `grid_extension.py` — computes and justifies the extended ladders.
- `add_recoverable_mass_multiband.py`, `calibrate_floors.py` — recoverable-mass columns and detection-floor calibration.
- `radmc3d/` — `userdef_module.f90`, `+radmc3d.inp`, opacity, wavelength, and external-source files.

### inject/
- `run_inject_v3.py` — injection driver; all settings are in the configuration files.
- `inject_config_v3.py`, `_flat.py`, `_iso.py`, `_flat_iso.py` — the four variants (as placed, flattened background, isolated, both).
- `inject.py` — placement primitives.
- `build_joint_target_v3.py`, `joint_target_Aquila_v3.txt` — target table of 1030 real Aquila sources.
- `dryrun_placement_v3.py` — selection and placement without writing images.
- `smoke_test_v3.py` — end-to-end test on synthetic inputs.

### tests/
- `run_validation.py` — all four validation tests (node leave-one-out, subdivision grid, background scan, population injection).
- `analyse_quad_injection.py`, `paired_injection_v3.py`, `analyse_injection_v3.py` — analysis of the four-variant injections, core by core.
- `check_injection_match_v3.py` — injected versus real population match.
- `sweep_sed_bands.py` — SED band-threshold sweep against known truth.
- `compare_size_measures.py` — model half-maximum width versus getsf AFWHM.
- `loo_4axis_full_v2.py` — standalone leave-one-node-out test.
- `run_clouds.py` — applies the correction to the seven clouds.
- `analyse_aquila_sim.py` — re-analysis of the Koenyves et al. (2015) simulation.

### hicorrect/
- `hicorrect.for` — reference implementation (distributed with getsf).
- `hicorrect.py` — Python version, kept identical to the Fortran.
