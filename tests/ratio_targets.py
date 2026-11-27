#!/usr/bin/env python3
"""
ratio_targets.py -- can physical quantities be recovered more accurately as
ratios to their observational counterparts than in absolute terms, and how does
the accuracy depend on the footprint ratio the observer actually measured?

TWO QUESTIONS
-------------
1. The mass correction succeeds because it predicts a ratio.  The corrector
   does not predict the mass, which the four observables determine only to
   0.114 decimal exponents, but the factor C_M = M_BE / M_ex by which the
   measured mass must be multiplied, which they determine to 0.027.  The
   measured mass carries most of the information and the prediction supplies
   only the correction to it.  The same construction is tried here for three
   further quantities:

     C_R     = R_BE / R_dec, the ratio of the true Bonnor-Ebert radius to the
               radius obtained by Gaussian deconvolution of the measured width
               on the surface-density band;
     C_T     = Td_avg / T_SED, the ratio of the mass-averaged dust temperature
               to the temperature returned by the fit to the spectral energy
               distribution;
     C_alpha = alpha_true / alpha_obs, the ratio of the true stability
               parameter to the one the observer forms.

   Here alpha = M_crit / M with M_crit = 2.4 R c_s^2 / G.  The observer's
   alpha_obs uses the deconvolved radius, the fitted spectral temperature and
   the measured mass; the true alpha_true uses the model's own radius,
   isothermal temperature and mass.  For a Bonnor-Ebert sphere alpha_true is
   below unity exactly when the sphere is super-critical, so the physically
   exact classification threshold on the corrected alpha is 1, not the
   conventional 2, which exists only to absorb the uncertainties.

   Each ratio is compared with the direct prediction of the quantity itself.

2. Every row of the recovery table is a measurement of a model through an
   aperture of a given size, and the footprint ratio phi that results is what
   an observer measures.  No representative value is imposed here: all samples
   are scored, and the classification accuracy is reported in bins of phi, so
   that the accuracy can be read at the values real catalogs exhibit, which lie
   between about 1.83 and 2.08.

METHOD
------
Local linear regression on the four validated observables, with the node under
test withheld, exactly as for the mass correction.

USAGE
    python3 ratio_targets.py
    python3 ratio_targets.py --targets xi_max C_R C_alpha
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


def loo(X, key, y):
    """Leave-one-node-out prediction of log10(y) at every sample."""
    pred = np.full(len(y), np.nan)
    for k in np.unique(key):
        out = key == k
        keep = ~out & np.isfinite(y)
        if keep.sum() < 50:
            continue
        est = mc.LocalLinear(X[keep], lg(y[keep]))
        pred[out] = est(X[out])
    return 10 ** pred


def summarize(label, pred, true):
    g = np.isfinite(pred) & np.isfinite(true) & (pred > 0) & (true > 0)
    lr = lg(pred[g] / true[g])
    print('  %-46s n=%6d  typical error %6.3f dex  scatter %5.2f  '
          'within 1.5x %3.0f%%'
          % (label, int(g.sum()), np.median(np.abs(lr)), 10 ** lr.std(),
             100 * np.mean((pred[g] / true[g] > 1 / 1.5)
                           & (pred[g] / true[g] < 1.5))))


def classify(label, truth_super, pred_super, phi):
    acc = 100 * np.mean(truth_super == pred_super)
    sens = 100 * pred_super[truth_super].mean()
    spec = 100 * (~pred_super[~truth_super]).mean()
    line = ('  %-42s  overall %5.1f%%  found %5.1f%%  kept %5.1f%%  |'
            % (label, acc, sens, spec))
    for lo, hi in zip(PHI_BINS[:-1], PHI_BINS[1:]):
        m = (phi >= lo) & (phi < hi)
        line += ' %5.1f' % (100 * np.mean(truth_super[m] == pred_super[m])
                            if m.sum() > 20 else np.nan)
    print(line)


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
    phi = S['cfoot'][ok]
    mrec, cmtrue = S['mrec'][ok], S['cm'][ok]
    node = np.asarray(S['node'])[ok]
    key = node[:, 0] * 100000 + node[:, 1] * 1000 + node[:, 2]

    Gr = read_grid(a.grid)
    ckeys, C = read_named(a.catalog, CAT_COLS)
    cat = {k: {n: C[n][i] for n in CAT_COLS} for i, k in enumerate(ckeys)}
    get = lambda src, n: np.array([src[k][n] if k in src else np.nan
                                   for k in key])                  # noqa: E731

    xi = get(Gr, 'xi_max')
    mbe = get(Gr, 'M_BE')
    tbe = get(Gr, 'T_BE')
    rbe_as = get(cat, 'R_BE_as')
    tdavg = get(cat, 'Td_avg')
    tsed = get(cat, 'T_SED4bsl')
    t3 = get(cat, 'T_SED3bsl')
    tsed = np.where(np.isfinite(tsed) & (tsed > 0), tsed, t3)
    fw = get(cat, 'FWHMSDbsl')

    dec_as = np.sqrt(np.clip(fw ** 2 - BEAM_AS ** 2, 0, None))
    scale = (mc.GRID_DISTANCE_PC / 206264.806) * PC
    a_true = (2.4 * (rbe_as * scale) * (KB * tbe / (MU * MH))
              / G_CGS / MSUN) / mbe
    a_obs = (2.4 * (0.5 * dec_as * scale) * (KB * tsed / (MU * MH))
             / G_CGS / MSUN) / mrec

    good = (np.all(np.isfinite(X), axis=1) & np.isfinite(xi)
            & np.isfinite(a_true) & (a_true > 0) & np.isfinite(a_obs)
            & (a_obs > 0))
    X, phi, key, xi = X[good], phi[good], key[good], xi[good]
    mrec, cmtrue, mbe = mrec[good], cmtrue[good], mbe[good]
    rbe_as, dec_as = rbe_as[good], dec_as[good]
    tdavg, tsed, a_true, a_obs = (tdavg[good], tsed[good], a_true[good],
                                  a_obs[good])
    print('samples scored: %d on %d nodes\n' % (len(key), len(np.unique(key))))

    print('Accuracy of the inversion, every sample scored, node withheld.\n'
          'The upper block predicts quantities directly; the lower block '
          'predicts the ratio\nof the quantity to its observational '
          'counterpart, which is then applied to the\nmeasured value.\n')
    for lab, y, truth in (('truncation radius xi_max', xi, xi),
                          ('Bonnor-Ebert radius R_BE', rbe_as, rbe_as),
                          ('mass-averaged temperature Td_avg', tdavg, tdavg),
                          ('stability parameter alpha', a_true, a_true)):
        summarize(lab, loo(X, key, y), truth)
    print()
    for lab, ratio, applied, truth in (
            ('mass factor C_M = M_BE/M_ex', cmtrue, mrec, mbe),
            ('size factor C_R = R_BE/R_dec', rbe_as / (0.5 * dec_as),
             0.5 * dec_as, rbe_as),
            ('temperature factor C_T = Td_avg/T_SED', tdavg / tsed, tsed,
             tdavg),
            ('stability factor C_alpha = a_true/a_obs', a_true / a_obs, a_obs,
             a_true)):
        pr = loo(X, key, ratio)
        summarize(lab, pr * applied, truth)

    print('\nStability classification, every sample scored.  "found" is the '
          'percentage of\ntruly super-critical samples so called, "kept" the '
          'percentage of truly\nsub-critical samples so called.  The columns '
          'after the bar give the accuracy in\nbins of the measured footprint '
          'ratio phi: %s.\n'
          % ', '.join('%.1f-%.1f' % (l, h)
                      for l, h in zip(PHI_BINS[:-1], PHI_BINS[1:])))
    ts = xi > XI_CRIT
    classify('conventional, alpha_obs <= 2', ts, a_obs <= 2.0, phi)
    classify('conventional, alpha_obs <= 1', ts, a_obs <= 1.0, phi)
    classify('inversion of xi_max', ts, loo(X, key, xi) > XI_CRIT, phi)
    ca = loo(X, key, a_true / a_obs)
    classify('corrected alpha = C_alpha x alpha_obs <= 1', ts,
             (ca * a_obs) <= 1.0, phi)
    print('\n  fraction of samples truly super-critical: %.1f%%'
          % (100 * ts.mean()))
    print('  median measured footprint ratio phi: %.2f' % np.median(phi))


if __name__ == '__main__':
    main()
