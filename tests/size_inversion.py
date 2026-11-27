"""Can the model grid serve as a deconvolution, that is, recover the true
truncation radius of a source from what an observer measures?

The question is not answered by predicting R_BE from the four observables used
by the mass corrector, because on the model side the fourth of those IS the
angle subtended by the truncation radius, so such a test is circular.  Here the
predictor is given only quantities an observer measures:

    Sigma        background surface density on the surface-density band,
    zeta         concentration,
    phi          footprint factor,
    r_as         the footprint radius itself, in arcseconds,
    w_sd         the measured width of the profile on the surface-density band,

with r_as reconstructed per sample as phi w_sd / 2, and is asked for the true
truncation radius R_BE in arcseconds.  One node is withheld at a time.

The grid-based estimate is compared with the three estimators an observer might
otherwise use: the footprint radius itself, half the measured width, and half
the Gaussian-deconvolved width, the last being standard practice.
"""
import sys
sys.path.insert(0, 'python')
import numpy as np
import mass_correction_pipeline_6 as mc
lg = np.log10
BEAM = 13.5


def q(label, ratio):
    g = np.isfinite(ratio) & (ratio > 0)
    lr = lg(ratio[g])
    print('  %-44s n=%5d  median %6.3f  scatter %5.2f  within 1.5x %3.0f%%  '
          'typical error %6.3f dex'
          % (label, int(g.sum()), np.median(ratio[g]), 10 ** lr.std(),
             100 * np.mean((ratio[g] > 1 / 1.5) & (ratio[g] < 1.5)),
             np.median(np.abs(lr))), flush=True)


mc.VERBOSE = 0
S = mc.build_recovery_samples('cats/bes_model_grid_final3_catalog',
                              'cats/bes_model_grid_final3_recovery_tables')
ok = np.isfinite(S['rbe_pc']) & (S['rbe_pc'] > 0)
rbe = S['rbe_pc'][ok] * 206264.806 / mc.GRID_DISTANCE_PC     # arcsec, truth
w = S['w_sd'][ok]
phi = S['cfoot'][ok]
r_as = phi * w / 2.0
X = np.column_stack([lg(S['sig'][ok]), lg(S['conc'][ok]), lg(phi), lg(r_as),
                     lg(w)])
node = np.asarray(S['node'])[ok]
key = node[:, 0] * 100000 + node[:, 1] * 1000 + node[:, 2]
g = np.all(np.isfinite(X), axis=1) & np.isfinite(rbe) & (rbe > 0)
X, key, rbe, w, r_as = X[g], key[g], rbe[g], w[g], r_as[g]
print('%d samples on %d nodes\n' % (len(key), len(np.unique(key))), flush=True)

pred = np.full(len(key), np.nan)
for k in np.unique(key):
    o = key == k
    pred[o] = mc.LocalLinear(X[~o], lg(rbe[~o]))(X[o])
pred = 10 ** pred

dec = np.sqrt(np.clip(w ** 2 - BEAM ** 2, 0, None))
print('Ratio of each estimate of the truncation radius to its true value:\n')
q('half the Gaussian-deconvolved width', 0.5 * dec / rbe)
q('half the measured width', 0.5 * w / rbe)
q('the footprint radius', r_as / rbe)
q('grid inversion, node withheld', pred / rbe)
