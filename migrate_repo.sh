#!/bin/bash
# migrate_repo.sh -- copy the files listed in MIGRATION.md from the old repository
# amenshch/Aquila-cutout-test into the local copy of the new repository, then commit
# and push them.
#
# Usage:
#   ./migrate_repo.sh <new-repo-folder>            copy, commit, push
#   ./migrate_repo.sh <new-repo-folder> --no-push  copy and commit only
#
# The old repository is downloaded into a temporary folder; only the listed files are
# fetched, not the whole repository.  Existing files in the new repository with the
# same names are overwritten.  Works with the bash 3.2 that comes with macOS.

set -e
OLD_URL="https://github.com/amenshch/Aquila-cutout-test.git"

NEW="$1"
PUSH=yes
[ "$2" = "--no-push" ] && PUSH=no
if [ -z "$NEW" ] || [ ! -d "$NEW/.git" ]; then
  echo "Give the folder of the local copy of the new repository (it must contain .git)."
  exit 1
fi
NEW=$(cd "$NEW" && pwd)

# Mapping: old path (file, folder, or pattern) | new folder.
# A trailing "/" in the new path means "copy into this folder".
MAP='
python/+Method/mass_correction_pipeline_6.py|method/
python/+Method/fit_mass.py|method/
python/+Method/getsf_columns.py|method/
python/+Method/getsf_bg.py|method/
python/+Method/read_getsf_sed.py|method/
python/+Grid/make_bes_grid3.py|grid/
python/+Grid/make_bes_grid3_subdivision.py|grid/
python/+Grid/tabulate_recovery.py|grid/
python/+Grid/flag_hole_models.py|grid/
python/+Grid/grid_extension.py|grid/
python/+Grid/add_recoverable_mass_multiband.py|grid/
python/+Grid/calibrate_floors.py|grid/
hires/radmc3d/*|grid/radmc3d/
python/+Inject/run_inject_v3.py|inject/
python/+Inject/inject.py|inject/
python/+Inject/inject_config_v3.py|inject/
python/+Inject/inject_config_v3_flat.py|inject/
python/+Inject/inject_config_v3_iso.py|inject/
python/+Inject/inject_config_v3_flat_iso.py|inject/
python/+Inject/build_joint_target_v3.py|inject/
python/+Inject/joint_target_Aquila_v3.txt|inject/
python/+Inject/dryrun_placement_v3.py|inject/
python/+Inject/smoke_test_v3.py|inject/
python/+Tests/run_validation.py|tests/
python/+Tests/analyse_quad_injection.py|tests/
python/+Tests/analyse_injection_v3.py|tests/
python/+Tests/paired_injection_v3.py|tests/
python/+Tests/check_injection_match_v3.py|tests/
python/+Tests/sweep_sed_bands.py|tests/
python/+Tests/compare_size_measures.py|tests/
python/+Tests/loo_4axis_full_v2.py|tests/
python/+Tests/loo4d_full_results_v2.csv|tests/
python/+Tests/run_clouds.py|tests/
python/+Tests/analyse_aquila_sim.py|tests/
python/+Figs/aa_plotstyle.py|figures/
python/+Figs/dkbluered.py|figures/
python/+Figs/make_flowchart_mpl.py|figures/
python/+Figs/make_fig_grid_coverage3.py|figures/
python/+Figs/make_fig_degeneracy.py|figures/
python/+Figs/make_fig_factor_decomposition.py|figures/
python/+Figs/make_fig_footprint_model.py|figures/
python/+Figs/make_fig_validation.py|figures/
python/+Figs/make_fig_capping_tradeoff.py|figures/
python/+Figs/fig_sim_vs_real.py|figures/
python/+Misc/cmf_mle.py|misc/
python/+Misc/check_separation.py|misc/
python/+Misc/estimator_compare.py|misc/
python/+Misc/make_compare_pseudoreplication.py|misc/
fortran/hicorrect.for|hicorrect/
python/+Hires/hicorrect.py|hicorrect/
cats/bes_grid_final3.txt|data/grid/
cats/bes_grid_final3_subdivision.txt|data/grid/
cats/bes_model_grid_final3_catalog|data/grid/
cats/bes_model_grid_final3_recovery_tables|data/grid/
cats/bes_model_grid_final3_subdivision_catalog|data/grid/
cats/bes_model_grid_final3_subdivision_recovery_tables|data/grid/
cats/bes_model_grid_final2_profiles|data/grid/
hires/profiles_final3.txt.gz|data/grid/
hires/bes_grid_G0=*.txt|data/grid_G0/
hires/harvest_calib3_G0=*.txt|data/grid_G0/
Aquila-Guoyin/*|data/clouds/Aquila/
California-Guoyin/*|data/clouds/California/
CygnusX-Guoyin/*|data/clouds/CygnusX/
Ophiuchus-Guoyin/*|data/clouds/Ophiuchus/
OrionA-Guoyin/*|data/clouds/OrionA/
Scorpius-Guoyin/*|data/clouds/Scorpius/
W3W4W5-Guoyin/*|data/clouds/W3W4W5/
Aquila/*|data/clouds/Aquila-Konyves2015/
Aquila-injections/260810_*|data/injections/
Aquila-injections/260811_*|data/injections/
Aquila-injections/260812_*|data/injections/
Aquila-injections/260813_*|data/injections/
Aquila-injections/260814_*|data/injections/
tests/bgscan_*|data/injections/bgscan/
Aquila-sim/*|data/injections/Aquila-sim/
latex/rt_mass_correction_method.tex|latex/
latex/rt_mass_correction.tex|latex/
latex/rt_mass_correction.bib|latex/
latex/aa.cls|latex/
latex/aa.bst|latex/
latex/lineno.sty|latex/
latex/linenoaa.sty|latex/
plots/fig_method_flowchart.pdf|latex/
plots/fig_grid_coverage.pdf|latex/
plots/fig_footprint_model.pdf|latex/
plots/fig_degeneracy.pdf|latex/
plots/fig_validation.pdf|latex/
plots/fig_factor_decomposition.pdf|latex/
plots/fig_capping_tradeoff.pdf|latex/
plots/fig_fit_validation.pdf|latex/
plots/fig_gridconv.pdf|latex/
plots/completeness_limit.pdf|latex/
docs/Summary_from_chat_*|docs/
docs/recovery_notes.md|docs/
docs/injection_test_notes.md|docs/
docs/MANIFEST_raw_injection.md|docs/
'

TMP=$(mktemp -d /tmp/oldrepo.XXXXXX)
echo "Downloading the file list of the old repository ..."
git clone --quiet --depth 1 --filter=blob:none --no-checkout "$OLD_URL" "$TMP/old"
cd "$TMP/old"

# Fetch only the listed paths.
PATHS=$(echo "$MAP" | sed '/^$/d' | cut -d'|' -f1)
echo "Fetching the listed files (this can take a few minutes for the large tables) ..."
echo "$PATHS" | tr '\n' '\0' | xargs -0 git checkout HEAD -- 2>/dev/null || true

# Copy.
N=0; MISSING=""
echo "$MAP" | sed '/^$/d' > "$TMP/map"
while IFS='|' read -r SRC DST; do
  mkdir -p "$NEW/$DST"
  FOUND=no
  for F in $SRC; do
    if [ -e "$F" ]; then
      cp -R "$F" "$NEW/$DST"
      FOUND=yes; N=$((N+1))
    fi
  done
  [ "$FOUND" = no ] && MISSING="$MISSING $SRC"
done < "$TMP/map"

cd "$NEW"
# Ignore Mac and LaTeX by-products.
if [ ! -f .gitignore ]; then
  printf '.DS_Store\n*.aux\n*.log\n*.blg\n*.bbl\n*.synctex.gz\n*.out\n__pycache__/\n' > .gitignore
fi

echo "Copied $N files or folders into $NEW"
[ -n "$MISSING" ] && echo "NOT FOUND in the old repository:$MISSING"

git add -A
git commit --quiet -m "Files migrated from Aquila-cutout-test (see MIGRATION.md)" && echo "Committed."
if [ "$PUSH" = yes ]; then
  git push --quiet && echo "Pushed to GitHub."
fi
rm -rf "$TMP"
