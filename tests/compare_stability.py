#!/usr/bin/env python3
"""
compare_stability.py -- how well can a core be classified as gravitationally
sub-critical or super-critical, and how well can its mass-averaged temperature
be recovered, from quantities an observer actually measures?

THE COMPARISON
--------------
Two methods are scored on the same model grid, where the answer is known.

  (a) The conventional criterion.  The observer forms the critical
      Bonnor-Ebert mass M_crit = 2.4 R c_s^2 / G from a radius R obtained by
      Gaussian deconvolution of the measured full width at half maximum on the
      surface-density band, and a sound speed evaluated at the temperature
      returned by the fit to the spectral energy distribution.  The core counts
      as bound when alpha_BE = M_crit / M does not exceed 2.  This replicates
      the observational procedure exactly, including its known weakness: simple
      Gaussian deconvolution of a centrally condensed profile is inaccurate by
      large factors (Men'shchikov 2023).

  (b) The grid inversion.  The four observables already validated for the mass
      correction -- the background surface density, the concentration, the
      footprint ratio, and the angle subtended by the truncation radius -- are
      inverted onto the dimensionless truncation radius xi_max by local linear
      regression, with the node under test withheld, and the core counts as
      super-critical when the inferred xi_max exceeds 6.451.

Both are scored against the model's true stability, which is set by its own
xi_max, on the same nodes and at the same measurement point, namely the
truncation radius at which the footprint ratio is closest to the grid's median
value of 1.84.  Method (a) is scored twice, once on the reported mass and once
on the corrected mass, to separate the effect of the mass correction from the
effect of the size and temperature estimates.

THE TEMPERATURE
---------------
The same inversion is applied to two temperatures: T_BE, the isothermal
temperature of the Bonnor-Ebert model, and Td_avg, the mass-averaged dust
temperature of the model as computed by the radiative transfer.  The second is
the physically meaningful one and is not directly observable: the temperature
an observer obtains from the spectral energy distribution is weighted by
emission and is biased warm.  As a reference, the error of the fitted spectral
temperature against the mass-averaged temperature is reported alongside.

USAGE
    python3 compare_stability.py
    python3 compare_stability.py --rep-foot 1.84

Inputs expected in the current directory:
    python/mass_correction_pipeline_6.py, python/getsf_columns.py
    cats/bes_model_grid_final3_catalog
    cats/bes_model_grid_final3_recovery_tables
    cats/bes_grid_final3.txt
"""
import argparse
import sys

sys.path.insert(0, 'python')
import numpy as np                                     # noqa: E402
import mass_correction_pipeline_6 as mc                # noqa: E402

lg = np.log10
XI_CRIT = 6.451
BEAM_AS = 13.5           # beam of the surface-density band, arcsec
G_CGS = 6.674e-8
KB = 1.3807e-16
MH = 1.6726e-24
MU = 2.33
MSUN = 1.989e33
PC = 3.086e18
ALPHA_BOUND = 2.0        # alpha_BE at or below which a core counts as bound

# Column positions, 1-based as printed in the catalog header, of the quantities
# taken from the 67-column model catalog.
CAT_COLS = {'T_BE': 6, 'M_BE': 8, 'R_BE_as': 9, 'Td_avg': 12,
            'T_SED4bsl': 17, 'T_SED3bsl': 23, 'FWHMSDbsl': 61}


def read_named(path, cols):
    """Read named columns from a whitespace catalog with '#' comment lines."""
    out = {n: [] for n in cols}
    keys = []
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


def score(name, truth_super, pred_super):
    n = len(truth_super)
    acc = 100 * np.mean(truth_super == pred_super)
    sens = 100 * pred_super[truth_super].mean() if truth_super.any() else np.nan
    spec = (100 * (~pred_super[~truth_super]).mean()
            if (~truth_super).any() else np.nan)
    print('  %-42s n=%5d  accuracy %5.1f%%  super-critical found %5.1f%%  '
          'sub-critical kept %5.1f%%' % (name, n, acc, sens, spec))
    return acc


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--catalog', default='cats/bes_model_grid_final3_catalog')
    ap.add_argument('--rectab',
                    default='cats/bes_model_grid_final3_recovery_tables')
    ap.add_argument('--grid', default='cats/bes_grid_final3.txt')
    ap.add_argument('--rep-foot', type=float, default=1.84)
    a = ap.parse_args()
    mc.VERBOSE = 0

    S = mc.build_recovery_samples(a.catalog, a.rectab)
    ok = np.isfinite(S['rbe_pc']) & (S['rbe_pc'] > 0)
    theta = S['rbe_pc'][ok] * 206264.806 / mc.GRID_DISTANCE_PC
    X = np.column_stack([lg(S['sig'][ok]), lg(S['conc'][ok]),
                         lg(S['cfoot'][ok]), lg(theta)])
    cfoot = S['cfoot'][ok]
    mrec = S['mrec'][ok]
    cm = S['cm'][ok]
    node = np.asarray(S['node'])[ok]
    key = node[:, 0] * 100000 + node[:, 1] * 1000 + node[:, 2]

    Gr = read_grid(a.grid)
    ckeys, C = read_named(a.catalog, CAT_COLS)
    cat = {k: {n: C[n][i] for n in CAT_COLS} for i, k in enumerate(ckeys)}

    xi = np.array([Gr[k]['xi_max'] if k in Gr else np.nan for k in key])
    good = np.isfinite(xi) & np.all(np.isfinite(X), axis=1)
    X, cfoot, mrec, cm, key, xi = (X[good], cfoot[good], mrec[good],
                                   cm[good], key[good], xi[good])

    # One representative measurement per node.
    uniq = np.unique(key)
    rep = np.array([np.where(key == k)[0][
        np.argmin(np.abs(cfoot[np.where(key == k)[0]] - a.rep_foot))]
        for k in uniq])
    print('nodes = %d, samples = %d, representative footprint ratio = %.2f\n'
          % (len(uniq), len(key), a.rep_foot))

    # ---- inversion of xi_max, node withheld -------------------------------
    pred_xi = np.full(len(uniq), np.nan)
    for i, k in enumerate(uniq):
        keep = key != k
        est = mc.LocalLinear(X[keep], lg(xi[keep]))
        pred_xi[i] = float(est(X[rep[i]:rep[i] + 1])[0])
    pred_xi = 10 ** pred_xi

    # ---- conventional criterion -------------------------------------------
    fw = np.array([cat[k]['FWHMSDbsl'] if k in cat else np.nan for k in uniq])
    Tsed = np.array([cat[k]['T_SED4bsl'] if k in cat else np.nan for k in uniq])
    T3 = np.array([cat[k]['T_SED3bsl'] if k in cat else np.nan for k in uniq])
    Tsed = np.where(np.isfinite(Tsed) & (Tsed > 0), Tsed, T3)
    Tdavg = np.array([cat[k]['Td_avg'] if k in cat else np.nan for k in uniq])
    m_ex = mrec[rep]
    m_corr = mrec[rep] * cm[rep]

    dec = np.sqrt(np.clip(fw ** 2 - BEAM_AS ** 2, 0, None))
    R = 0.5 * dec * (mc.GRID_DISTANCE_PC / 206264.806) * PC
    Mcrit = 2.4 * R * (KB * Tsed / (MU * MH)) / G_CGS / MSUN

    truth_super = xi[rep] > XI_CRIT
    val = (np.isfinite(Mcrit) & (Mcrit > 0) & np.isfinite(m_ex) & (m_ex > 0)
           & np.isfinite(Tsed) & (Tsed > 0) & np.isfinite(pred_xi))
    print('Stability classification, scored against the true xi_max of each '
          'model.\n"super-critical found" is the percentage of truly '
          'super-critical nodes so called;\n"sub-critical kept" is the '
          'percentage of truly sub-critical nodes so called.\n')
    base = max(truth_super[val].mean(), 1 - truth_super[val].mean())
    print('  truly super-critical: %.1f%% of the %d nodes scored, so always '
          'answering with the\n  commoner class would score %.1f%%\n'
          % (100 * truth_super[val].mean(), int(val.sum()), 100 * base))
    score('conventional criterion, reported mass',
          truth_super[val], (Mcrit / m_ex)[val] <= ALPHA_BOUND)
    score('conventional criterion, corrected mass',
          truth_super[val], (Mcrit / m_corr)[val] <= ALPHA_BOUND)
    score('grid inversion of xi_max, node withheld',
          truth_super[val], pred_xi[val] > XI_CRIT)

    dfw = dec / (2 * np.array([cat[k]['R_BE_as'] if k in cat else np.nan
                               for k in uniq]))
    v2 = np.isfinite(dfw) & (dfw > 0)
    print('\n  For reference, the ratio of the Gaussian-deconvolved diameter '
          'to the true\n  Bonnor-Ebert diameter has median %.2f and spans '
          '%.2f to %.2f between its\n  tenth and ninetieth percentiles over '
          '%d nodes.'
          % (np.median(dfw[v2]), *np.percentile(dfw[v2], [10, 90]),
             int(v2.sum())))

    # ---- temperature inversion --------------------------------------------
    print('\nInference of temperature from the same four observables, node '
          'withheld.\nThe error is the median absolute difference between the '
          'inferred and the true\nvalue, in kelvin, and as the median '
          'absolute base-ten logarithm of the ratio.\n')
    for label, tv in (('isothermal model temperature T_BE',
                       np.array([Gr[k]['T_BE'] for k in uniq])),
                      ('mass-averaged dust temperature Td_avg', Tdavg)):
        tsample = np.array([Gr[k]['T_BE'] if 'T' in label and
                            label.startswith('isothermal') else np.nan
                            for k in key])
        if label.startswith('isothermal'):
            ysamp = np.array([Gr[k]['T_BE'] for k in key])
        else:
            ysamp = np.array([cat[k]['Td_avg'] if k in cat else np.nan
                              for k in key])
        gg = np.isfinite(ysamp) & (ysamp > 0)
        pr = np.full(len(uniq), np.nan)
        for i, k in enumerate(uniq):
            keep = gg & (key != k)
            est = mc.LocalLinear(X[keep], lg(ysamp[keep]))
            pr[i] = float(est(X[rep[i]:rep[i] + 1])[0])
        pr = 10 ** pr
        m = np.isfinite(pr) & np.isfinite(tv) & (tv > 0)
        print('  %-40s n=%5d  median error %5.2f K  %6.3f dex'
              % (label, int(m.sum()), np.median(np.abs(pr[m] - tv[m])),
                 np.median(np.abs(lg(pr[m] / tv[m])))))
        del tsample
    m = np.isfinite(Tsed) & np.isfinite(Tdavg) & (Tdavg > 0) & (Tsed > 0)
    print('  %-40s n=%5d  median error %5.2f K  %6.3f dex'
          % ('fitted spectral temperature, as reference', int(m.sum()),
             np.median(np.abs(Tsed[m] - Tdavg[m])),
             np.median(np.abs(lg(Tsed[m] / Tdavg[m])))))
    print('  (the last line is the error an observer incurs by adopting the '
          'temperature\n  from the fit to the spectral energy distribution as '
          'the mass-averaged one)')


if __name__ == '__main__':
    main()
