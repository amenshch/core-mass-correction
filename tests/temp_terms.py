"""The temperature term of the mass correction, measured rather than inferred.

The pipeline reports a decomposition C_M = F_flux x F_bg x F_temp in which
F_temp is defined as the residual C_M/(F_flux F_bg), each factor being
separately interpolated, and is then clamped from below at unity.  It is
therefore not a measurement of any temperature effect.  Here the terms are
formed directly from quantities the model catalogue knows exactly:

  M_BE / M_SED4bsl   the ratio of the true mass to the mass an observer would
                     derive from a modified-blackbody fit to the whole model
                     with the true crater background removed: no truncation
                     loss and no background error remain, so this is the pure
                     bias of the spectral-energy-distribution mass, which is
                     where any temperature-structure effect must appear;
  M_BE / M_SED4bs    the same with the interpolated background, so the
                     difference between the two is the background term;
  T_SED / Td_avg     the corresponding temperature ratio.

The residual C_M/(F_flux F_bg) is also formed from the exact per-sample values,
to show what it does and does not contain.
"""
import sys
sys.path.insert(0, 'python')
import numpy as np
import mass_correction_pipeline_6 as mc
lg = np.log10
COLS = {'M_BE': 8, 'Td_avg': 12, 'T4bsl': 17, 'M4bsl': 18, 'T4bs': 14,
        'M4bs': 15, 'T3bsl': 23, 'M3bsl': 24, 'T3bs': 20, 'M3bs': 21}


def read(path, cols):
    keys, out = [], {n: [] for n in cols}
    for l in open(path):
        if l.startswith('#') or not l.strip():
            continue
        f = l.split(); keys.append((int(f[1]), int(f[2]), int(f[3])))
        for n, i in cols.items():
            try: out[n].append(float(f[i - 1]))
            except Exception: out[n].append(np.nan)
    return keys, {n: np.array(v) for n, v in out.items()}


def q(label, v):
    g = np.isfinite(v) & (v > 0)
    if g.sum() < 10:
        print('  %-46s too few' % label); return
    print('  %-46s n=%5d  median %6.3f  tenth %6.3f  ninetieth %6.3f'
          % (label, int(g.sum()), np.median(v[g]),
             *np.percentile(v[g], [10, 90])))


k, C = read('cats/bes_model_grid_final3_catalog', COLS)
print('Per model, over the whole grid catalogue:\n')
q('true mass / SED mass, 4 bands, crater background',
  C['M_BE'] / C['M4bsl'])
q('true mass / SED mass, 4 bands, interpolated bg', C['M_BE'] / C['M4bs'])
q('true mass / SED mass, 3 bands, crater background',
  C['M_BE'] / C['M3bsl'])
q('SED temperature / mass-averaged temperature, 4 bands, crater',
  C['T4bsl'] / C['Td_avg'])
q('SED temperature / mass-averaged temperature, 4 bands, interp',
  C['T4bs'] / C['Td_avg'])

mc.VERBOSE = 0
S = mc.build_recovery_samples('cats/bes_model_grid_final3_catalog',
                              'cats/bes_model_grid_final3_recovery_tables')
print('\nPer training sample, from the exact per-sample values:\n')
q('C_M = M_BE / M_reported', S['cm'])
q('F_flux = M_SED3bs / M_reported', S['f_flux'])
q('F_bg = FXSDbsl / FXSDbs', S['f_bg'])
q('residual C_M / (F_flux F_bg), unclamped',
  S['cm'] / (S['f_flux'] * S['f_bg']))
r = S['cm'] / (S['f_flux'] * S['f_bg'])
g = np.isfinite(r) & (r > 0)
print('  fraction of samples with the residual below unity, which the '
      'pipeline clamps to 1: %.1f%%' % (100 * np.mean(r[g] < 1)))
