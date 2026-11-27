"""
paths.py -- locations of all data files of the repository, in one place.

Every script reads its default input files from here, so that a change of grid
or a moved folder needs one edit only.  All paths are absolute, built from the
location of this file, so scripts work from any current directory.

Usage in a script located in any folder of the repository:

    import os, sys
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                    '..', 'method'))
    import paths

The model grid in use is final3 (2443 models, Summary of chat 13); its
subdivision grid lies at the geometric midpoints of the final3 ladders and is
never used for training.
"""
import os

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# ---------------------------------------------------------------- model grid
GRID_DIR = os.path.join(REPO, 'data', 'grid')
GRID_NODES = os.path.join(GRID_DIR, 'bes_grid_final3.txt')
GRID_CATALOG = os.path.join(GRID_DIR, 'bes_model_grid_final3_catalog')
GRID_RECTAB = os.path.join(GRID_DIR, 'bes_model_grid_final3_recovery_tables')
GRID_PROFILES = os.path.join(GRID_DIR, 'profiles_final3.txt.gz')

SUBDIV_NODES = os.path.join(GRID_DIR, 'bes_grid_final3_subdivision.txt')
SUBDIV_CATALOG = os.path.join(GRID_DIR, 'bes_model_grid_final3_subdivision_catalog')
SUBDIV_RECTAB = os.path.join(GRID_DIR,
                             'bes_model_grid_final3_subdivision_recovery_tables')

# The recoverable-mass catalog read only by the legacy "loo" mode of the
# pipeline exists for the final2 grid alone and is not in this repository.
GRID_CATREC = os.path.join(GRID_DIR, 'bes_model_grid_final2_catalog_rec')

# Radial profiles of final2 nodes, still read by the injection code.
GRID_PROFILES_FINAL2 = os.path.join(GRID_DIR, 'bes_model_grid_final2_profiles')

# Grid ladder over the strength G_0 of the interstellar radiation field.
GRID_G0_DIR = os.path.join(REPO, 'data', 'grid_G0')

# -------------------------------------------------------------- observations
# One folder per cloud, holding the getsf catalogs (.cat, .add.cat, and the
# concatenated SED catalog).  "{cloud}" is replaced by the cloud name.
CLOUD_DIR = os.path.join(REPO, 'data', 'clouds', '{cloud}')
KONYVES2015_DIR = os.path.join(REPO, 'data', 'clouds', 'Aquila-Konyves2015')

# Clouds and adopted distances (pc), as in Zhang & Men'shchikov (2026).
CLOUDS = [('Scorpius', 130), ('Ophiuchus', 144), ('Aquila', 260),
          ('OrionA', 432), ('California', 470), ('CygnusX', 1150),
          ('W3W4W5', 1700)]

# ---------------------------------------------------------------- injections
INJECTION_DIR = os.path.join(REPO, 'data', 'injections')
BGSCAN_DIR = os.path.join(INJECTION_DIR, 'bgscan')
AQUILA_SIM_DIR = os.path.join(INJECTION_DIR, 'Aquila-sim')
TARGET_AQUILA = os.path.join(REPO, 'inject', 'joint_target_Aquila_v3.txt')

# ------------------------------------------------------------------- results
RESULTS_DIR = os.path.join(REPO, 'data', 'results')
CHAT14_RESULTS_DIR = os.path.join(REPO, 'tests', 'results_chat14')


def cloud_dir(cloud):
    """Folder of the getsf catalogs of one cloud."""
    return CLOUD_DIR.format(cloud=cloud)


if __name__ == '__main__':
    # Report which of the files named above exist.
    for name, value in sorted(globals().items()):
        if name.isupper() and isinstance(value, str) and '{' not in value:
            print('%-22s %-7s %s' % (name, 'ok' if os.path.exists(value) else 'MISSING',
                                     os.path.relpath(value, REPO)))
