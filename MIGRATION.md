# Migration from Aquila-cutout-test to the new mass-correction repository

This list maps files of the old repository (`amenshch/Aquila-cutout-test`, branch `main`, 2953 files) to the new layout.
Only files needed to continue the mass correction method are copied.
The old repository stays unchanged as an archive.
Paths are relative to the root of each repository.

Note on uploading through the GitHub web interface: a single file must be smaller than 25 MiB, and at most 100 files can be uploaded at once.
The two largest files below (`grid_final3_recovery_tables`, 24.1 MB, and `grid_final3_subdivision_recovery_tables`, 21.5 MB) are just below that limit.

## 1. Files to copy

### method/ — the correction pipeline

| Old path | New path |
|---|---|
| python/+Method/mass_correction_pipeline_6.py | method/mass_correction_pipeline_6.py |
| python/+Method/fit_mass.py | method/fit_mass.py |
| python/+Method/getsf_columns.py | method/getsf_columns.py |
| python/+Method/getsf_bg.py | method/getsf_bg.py |
| python/+Method/read_getsf_sed.py | method/read_getsf_sed.py |

### grid/ — model grid definition and recovery tables

| Old path | New path |
|---|---|
| python/+Grid/make_bes_grid3.py | grid/make_bes_grid3.py |
| python/+Grid/make_bes_grid3_subdivision.py | grid/make_bes_grid3_subdivision.py |
| python/+Grid/tabulate_recovery.py | grid/tabulate_recovery.py |
| python/+Grid/flag_hole_models.py | grid/flag_hole_models.py |
| python/+Grid/grid_extension.py | grid/grid_extension.py |
| python/+Grid/add_recoverable_mass_multiband.py | grid/add_recoverable_mass_multiband.py |
| python/+Grid/calibrate_floors.py | grid/calibrate_floors.py |
| hires/radmc3d/ (all 9 files) | grid/radmc3d/ |

### inject/ — injection experiments

| Old path | New path |
|---|---|
| python/+Inject/run_inject_v3.py | inject/run_inject_v3.py |
| python/+Inject/inject.py | inject/inject.py |
| python/+Inject/inject_config_v3.py | inject/inject_config_v3.py |
| python/+Inject/inject_config_v3_flat.py | inject/inject_config_v3_flat.py |
| python/+Inject/inject_config_v3_iso.py | inject/inject_config_v3_iso.py |
| python/+Inject/inject_config_v3_flat_iso.py | inject/inject_config_v3_flat_iso.py |
| python/+Inject/build_joint_target_v3.py | inject/build_joint_target_v3.py |
| python/+Inject/joint_target_Aquila_v3.txt | inject/joint_target_Aquila_v3.txt |
| python/+Inject/dryrun_placement_v3.py | inject/dryrun_placement_v3.py |
| python/+Inject/smoke_test_v3.py | inject/smoke_test_v3.py |

### tests/ — validation and analysis

| Old path | New path |
|---|---|
| python/+Tests/run_validation.py | tests/run_validation.py |
| python/+Tests/analyse_quad_injection.py | tests/analyse_quad_injection.py |
| python/+Tests/analyse_injection_v3.py | tests/analyse_injection_v3.py |
| python/+Tests/paired_injection_v3.py | tests/paired_injection_v3.py |
| python/+Tests/check_injection_match_v3.py | tests/check_injection_match_v3.py |
| python/+Tests/sweep_sed_bands.py | tests/sweep_sed_bands.py |
| python/+Tests/compare_size_measures.py | tests/compare_size_measures.py |
| python/+Tests/loo_4axis_full_v2.py | tests/loo_4axis_full_v2.py |
| python/+Tests/loo4d_full_results_v2.csv | tests/loo4d_full_results_v2.csv (its output; final2 grid) |
| python/+Tests/run_clouds.py | tests/run_clouds.py |
| python/+Tests/analyse_aquila_sim.py | tests/analyse_aquila_sim.py (a second copy exists in python/+Inject/; check which is current) |

### figures/ — paper figures

| Old path | New path |
|---|---|
| python/+Figs/aa_plotstyle.py | figures/aa_plotstyle.py |
| python/+Figs/dkbluered.py | figures/dkbluered.py |
| python/+Figs/make_flowchart_mpl.py | figures/make_flowchart_mpl.py |
| python/+Figs/make_fig_grid_coverage3.py | figures/make_fig_grid_coverage3.py |
| python/+Figs/make_fig_degeneracy.py | figures/make_fig_degeneracy.py |
| python/+Figs/make_fig_factor_decomposition.py | figures/make_fig_factor_decomposition.py |
| python/+Figs/make_fig_footprint_model.py | figures/make_fig_footprint_model.py |
| python/+Figs/make_fig_validation.py | figures/make_fig_validation.py |
| python/+Figs/make_fig_capping_tradeoff.py | figures/make_fig_capping_tradeoff.py |
| python/+Figs/fig_sim_vs_real.py | figures/fig_sim_vs_real.py |

### misc/

| Old path | New path |
|---|---|
| python/+Misc/cmf_mle.py | misc/cmf_mle.py |
| python/+Misc/check_separation.py | misc/check_separation.py |
| python/+Misc/estimator_compare.py | misc/estimator_compare.py |
| python/+Misc/make_compare_pseudoreplication.py | misc/make_compare_pseudoreplication.py |

### hicorrect/ — surface density correction (to be applied to grid and cloud images)

| Old path | New path |
|---|---|
| fortran/hicorrect.for | hicorrect/hicorrect.for (reference version) |
| python/+Hires/hicorrect.py | hicorrect/hicorrect.py |

### data/grid/

| Old path | New path |
|---|---|
| cats/bes_grid_final3.txt | data/grid/bes_grid_final3.txt |
| cats/bes_grid_final3_subdivision.txt | data/grid/bes_grid_final3_subdivision.txt |
| cats/bes_model_grid_final3_catalog | data/grid/bes_model_grid_final3_catalog |
| cats/bes_model_grid_final3_recovery_tables | data/grid/bes_model_grid_final3_recovery_tables |
| cats/bes_model_grid_final3_subdivision_catalog | data/grid/bes_model_grid_final3_subdivision_catalog |
| cats/bes_model_grid_final3_subdivision_recovery_tables | data/grid/bes_model_grid_final3_subdivision_recovery_tables |
| hires/profiles_final3.txt.gz | data/grid/profiles_final3.txt.gz |
| cats/bes_model_grid_final2_profiles | data/grid/bes_model_grid_final2_profiles (only if the injection code still reads it for the 95% flux radius) |

### data/grid_G0/ — radiation-field ladder (for the planned G_0 dependence)

| Old path | New path |
|---|---|
| hires/bes_grid_G0=*.txt (9 files) | data/grid_G0/ |
| hires/harvest_calib3_G0=*.txt (9 files) | data/grid_G0/ |

### data/clouds/

| Old path | New path |
|---|---|
| Aquila-Guoyin/, California-Guoyin/, CygnusX-Guoyin/, Ophiuchus-Guoyin/, OrionA-Guoyin/, Scorpius-Guoyin/, W3W4W5-Guoyin/ (3 files each) | data/clouds/<cloud>/ |
| Aquila/ (Koenyves et al. 2015 tables and ReadMe) | data/clouds/Aquila-Konyves2015/ |

### data/injections/

| Old path | New path |
|---|---|
| Aquila-injections/260810_* to 260814_* (18 runs; catalogs and truth tables) | data/injections/ |
| tests/bgscan_i* (7 background-scan runs) | data/injections/bgscan/ |
| Aquila-sim/ (Koenyves-era simulation, re-analyzed in the living document) | data/injections/Aquila-sim/ |

### latex/

| Old path | New path |
|---|---|
| latex/rt_mass_correction_method.tex | latex/rt_mass_correction_method.tex |
| latex/rt_mass_correction.tex | latex/rt_mass_correction.tex |
| latex/rt_mass_correction.bib | latex/rt_mass_correction.bib |
| latex/aa.cls, aa.bst, lineno.sty, linenoaa.sty | latex/ |
| plots/fig_method_flowchart.pdf, fig_grid_coverage.pdf, fig_footprint_model.pdf, fig_degeneracy.pdf, fig_validation.pdf, fig_factor_decomposition.pdf, fig_capping_tradeoff.pdf, fig_fit_validation.pdf, fig_gridconv.pdf, completeness_limit.pdf | latex/ (figures the two papers include) |

### docs/

| Old path | New path |
|---|---|
| docs/Summary_from_chat_1.txt to Summary_from_chat_9_addendum.md | docs/ |
| docs/recovery_notes.md, injection_test_notes.md, MANIFEST_raw_injection.md | docs/ |
| Project files Summary_from_chat_10.md, Summary_of_chat_12.txt, Summary_of_chat_13.md, Summary_of_chat_14.md | docs/ (no summary of chat 11 was made) |

## 2. Not migrated (remain in the archive)

- All hires development: `hires/` (except the items above), `python/+Hires/` (except hicorrect.py), the hires-correction drafts in `latex/26*`, and the FITS images in the repository root.
- Superseded hires-correction approaches: `python/+Method/hires_apply.py`, `hires_correct_levels.py`, `python/+Figs/make_letter_figures.py`, `validate_correction_models.py`.
- Superseded scripts: `python/old/`, `python/+Misc/patch_*.py`, duplicates in `plots/*.py`.
- Superseded grids: `cats/bes_grid_final.txt`, `bes_grid_final2.txt`, all `*final2*`, `*subdivision*` (except final3), `*params*` catalogs.
- Retired tests: the cSD injection files in `tests/` and the runs `Aquila-injections/260807_*`, `260808_*`.
- `getsf/+GETSF_WEBSITE`, `validation/Archive.tar.xz` (chat-12 archive).

## 3. Open points before migrating

1. The radiative-transfer model images of the final3 and subdivision grids are not in the repository (only about 20 nodes are under `hires/model_*`).
   Rebuilding the recovery tables from corrected surface density images needs all of them.
2. The injected-field images (the waveband and hires images of each injection run) are not in the repository either; only the catalogs and truth tables are.
   Rerunning getsf on corrected maps of injected fields needs them, or the injections must be regenerated.
3. `profile_fit_test.py` (the profile-fit fourth axis, cited as indicative in the method paper) was never preserved.
