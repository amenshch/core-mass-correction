#!/usr/bin/env python3
"""
branch_classifier.py -- two tests of how the stability of a core is best
decided from what an observer measures.

TEST 1: CLASSIFY THE BRANCH RATHER THAN REGRESS THE RADIUS
----------------------------------------------------------
The Bonnor-Ebert mass is a double-valued function of the dimensionless
truncation radius xi_max, taking the same value on either side of the critical
value 6.451.  Regressing the logarithm of xi_max therefore asks the estimator
to fit a quantity whose relation to the observables folds at the critical
point.  Here the estimator is given the branch itself as the target, the label
being 1 for super-critical and 0 for sub-critical, so that the local linear fit
returns a probability rather than a radius.  The same estimator, observables,
neighbour count and leave-one-node-out protocol are used, so any difference is
due to the change of target alone.

TEST 2: MEASURED SIZE INSTEAD OF DECONVOLVED SIZE
-------------------------------------------------
The conventional criterion forms the critical mass from a radius obtained by
Gaussian deconvolution of the measured width, which underestimates the true
Bonnor-Ebert radius by a median factor of about two and by up to a factor of
seven.  Since the critical mass scales linearly with the radius, the criterion
is correspondingly biased.  This test repeats the conventional criterion using
the measured half-width directly, with no deconvolution, to see whether the
customary threshold of 2 is in fact appropriate to an undeconvolved size.

Both tests score every sample, and report accuracy in bins of the measured
footprint ratio phi, since that is what an observer obtains and real catalogs
have median values between about 1.83 and 2.08.

USAGE
    python3 branch_classifier.py
"""
import argparse
import sys

sys.path.insert(0, 'python')
import numpy as np                                     # noqa: E402
import mass_correction_pipeline_6 as mc                # noqa: E402

lg = np.log10
XI_CRIT = 6.451
BEAM_AS = 13.5
G_CGS, KB, MH, MU = 6.674e-8, 1.3807e-16, 1.6726e-24, 2.33
MSUN, PC = 1.989e33, 3.086e18
CAT_COLS = {'T_BE': 6, 'M_BE': 8, 'R_BE_as': 9, 'Td_avg': 12,
            'T_SED4bsl': 17, 'T_SED3bsl': 23, 'FWHMSDbsl': 61}
PHI_BINS = [1.0, 1.4, 1.7, 1.9, 2.1, 2.5, 4.0]


def read_named(path, cols):
    keys, out = [], {n: [] for n in cols}
    for line in open(path):
        if line.startswith('#') or not line.strip():
            continue
        f = line.split()
        keys.append(int(f[1]) * 100000 + int(f[2]) * 1000 + int(f[3]))
        for n, c in cols.items():
            try:
                out[n].append(float(f[c - 1]))
            except (ValueError, IndexError):
                out[n].append(np.nan)
    return keys, {n: np.array(v) for n, v in out.items()}


def read_grid(path):
    rows = [l.split() for l in open(path) if l.strip()]
    col = {n: i for i, n in enumerate(rows[0])}
    out = {}
    for r in rows[1:]:
        key = (int(r[col['i']]) * 100000 + int(r[col['j']]) * 1000
               + int(r[col['k']]))
        d = {}
        for n, i in col.items():
            try:
                d[n] = float(r[i])
            except ValueError:
                d[n] = np.nan
        out[key] = d
    return out


def report(label, truth, pred, phi):
    line = ('  %-44s overall %5.1f%%  found %5.1f%%  kept %5.1f%%  |'
            % (label, 100 * np.mean(truth == pred),
               100 * pred[truth].mean(), 100 * (~pred[~truth]).mean()))
    for lo, hi in zip(PHI_BINS[:-1], PHI_BINS[1:]):
        m = (phi >= lo) & (phi < hi)
        line += ' %5.1f' % (100 * np.mean(truth[m] == pred[m])
                            if m.sum() > 20 else np.nan)
    print(line)


def best_threshold(truth, value, grid, phi=None, lo=None, hi=None):
    m = np.ones(len(truth), bool) if phi is None else ((phi >= lo)
                                                       & (phi < hi))
    b = max((100 * np.mean((value[m] <= t) == truth[m]), t) for t in grid)
    return b[1], b[0], int(m.sum())


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--catalog', default='cats/bes_model_grid_final3_catalog')
    ap.add_argument('--rectab',
                    default='cats/bes_model_grid_final3_recovery_tables')
    ap.add_argument('--grid', default='cats/bes_grid_final3.txt')
    a = ap.parse_args()
    mc.VERBOSE = 0

    S = mc.build_recovery_samples(a.catalog, a.rectab)
    ok = np.isfinite(S['rbe_pc']) & (S['rbe_pc'] > 0)
    theta = S['rbe_pc'][ok] * 206264.806 / mc.GRID_DISTANCE_PC
    X = np.column_stack([lg(S['sig'][ok]), lg(S['conc'][ok]),
                         lg(S['cfoot'][ok]), lg(theta)])
    phi, mrec = S['cfoot'][ok], S['mrec'][ok]
    node = np.asarray(S['node'])[ok]
    key = node[:, 0] * 100000 + node[:, 1] * 1000 + node[:, 2]

    Gr = read_grid(a.grid)
    ckeys, C = read_named(a.catalog, CAT_COLS)
    cat = {k: {n: C[n][i] for n in CAT_COLS} for i, k in enumerate(ckeys)}
    get = lambda s, n: np.array([s[k][n] if k in s else np.nan             # noqa
                                 for k in key])
    xi = get(Gr, 'xi_max')
    tsed = get(cat, 'T_SED4bsl')
    t3 = get(cat, 'T_SED3bsl')
    tsed = np.where(np.isfinite(tsed) & (tsed > 0), tsed, t3)
    fw = get(cat, 'FWHMSDbsl')

    good = (np.all(np.isfinite(X), axis=1) & np.isfinite(xi)
            & np.isfinite(fw) & (fw > 0) & np.isfinite(tsed) & (tsed > 0)
            & np.isfinite(mrec) & (mrec > 0))
    X, phi, key, xi, mrec, tsed, fw = (X[good], phi[good], key[good],
                                       xi[good], mrec[good], tsed[good],
                                       fw[good])
    truth = xi > XI_CRIT
    print('samples %d on %d nodes; truly super-critical %.1f%%; median '
          'footprint ratio %.2f\n'
          % (len(key), len(np.unique(key)), 100 * truth.mean(),
             np.median(phi)))

    scale = (mc.GRID_DISTANCE_PC / 206264.806) * PC
    cs2 = KB * tsed / (MU * MH)
    dec = np.sqrt(np.clip(fw ** 2 - BEAM_AS ** 2, 0, None))
    a_dec = 2.4 * (0.5 * dec * scale) * cs2 / G_CGS / MSUN / mrec
    a_mea = 2.4 * (0.5 * fw * scale) * cs2 / G_CGS / MSUN / mrec

    # ---- Test 1: branch classifier -----------------------------------------
    prob = np.full(len(key), np.nan)
    for k in np.unique(key):
        out = key == k
        est = mc.LocalLinear(X[~out], truth[~out].astype(float))
        prob[out] = est(X[out])
    ok2 = np.isfinite(prob)

    print('Stability classification, all samples, node withheld.  "found" is '
          'the percentage\nof truly super-critical samples so called, "kept" '
          'the percentage of truly\nsub-critical samples so called.  Columns '
          'after the bar: accuracy in bins of the\nmeasured footprint ratio '
          'phi = %s.\n'
          % ', '.join('%.1f-%.1f' % (l, h)
                      for l, h in zip(PHI_BINS[:-1], PHI_BINS[1:])))
    report('regression on xi_max (for reference)', truth[ok2],
           np.full(ok2.sum(), False), phi[ok2]) if False else None
    report('branch classifier, probability > 0.5', truth[ok2],
           prob[ok2] > 0.5, phi[ok2])
    t, acc, _ = best_threshold(~truth[ok2], -prob[ok2],
                               -np.arange(0.05, 0.96, 0.05))
    report('branch classifier, probability > %.2f' % -t, truth[ok2],
           prob[ok2] > -t, phi[ok2])

    # ---- Test 2: measured versus deconvolved size --------------------------
    grid = np.arange(0.2, 6.01, 0.05)
    report('conventional, deconvolved size, alpha <= 2', truth, a_dec <= 2,
           phi)
    report('conventional, measured size, alpha <= 2', truth, a_mea <= 2, phi)
    report('conventional, measured size, alpha <= 1', truth, a_mea <= 1, phi)
    for nm, av in (('deconvolved', a_dec), ('measured', a_mea)):
        t, acc, n = best_threshold(truth, av, grid)
        print('  best threshold with the %-12s size: alpha <= %.2f, '
              'accuracy %.1f%% over %d samples' % (nm, t, acc, n))
        for lo, hi in ((1.7, 1.9), (1.9, 2.1)):
            t2, acc2, n2 = best_threshold(truth, av, grid, phi, lo, hi)
            print('      within phi %.1f to %.1f: alpha <= %.2f, accuracy '
                  '%.1f%% over %d samples' % (lo, hi, t2, acc2, n2))
    print('\n  ratio of the measured half-width to the true Bonnor-Ebert '
          'radius: median %.2f'
          % np.median(0.5 * fw / get(cat, 'R_BE_as')[good]))
    print('  ratio of the deconvolved half-width to the same: median %.2f'
          % np.median(0.5 * dec / get(cat, 'R_BE_as')[good]))


if __name__ == '__main__':
    main()
