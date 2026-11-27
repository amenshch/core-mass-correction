#!/usr/bin/env python3
"""
extra_observables.py -- does adding the per-band sizes and the spectral colours
to the four validated observables improve the classification of a core as
sub-critical or super-critical?

RATIONALE
---------
The ratio of the critical Bonnor-Ebert mass to the measured mass cannot
separate the two branches, because the Bonnor-Ebert mass is double-valued in
the dimensionless truncation radius xi_max.  Only the shape of the source
breaks that degeneracy.  The four validated observables -- the background
surface density, the concentration, the footprint factor and the angle subtended
by the truncation radius -- carry shape information from the projected profile
alone, and with them the classification reaches 81-85 per cent at the footprint
factors real catalogues exhibit.  Two further kinds of information are present
in the models and have not been used:

  the per-band sizes, whose ratios measure how the apparent extent changes with
  wavelength, which is a direct expression of the internal temperature gradient
  and therefore of the concentration;

  the spectral colours, whose ratios measure the same gradient through the
  shape of the spectral energy distribution rather than through the images.

Both are dimensionless ratios and therefore distance-independent, as the four
validated observables are.  This script adds them, singly and together, and
scores the classification against the models' true stability, with the node
under test withheld throughout.  The comparison to beat is the accuracy in the
bins of the footprint factor that real observations reach.

CAVEAT
------
These additional observables have not been validated through the injection and
extraction chain, as the four have.  A gain measured here is a necessary but
not a sufficient condition for using them: it would have to be confirmed on
injected sources before it could be claimed.  A null result here, on the other
hand, settles the question without further work.

USAGE
    python3 extra_observables.py
"""
import argparse
import sys

sys.path.insert(0, 'python')
import numpy as np                                     # noqa: E402
import mass_correction_pipeline_6 as mc                # noqa: E402

lg = np.log10
XI_CRIT = 6.451
PHI_BINS = [1.0, 1.4, 1.7, 1.9, 2.1, 2.5, 4.0]
CAT_COLS = {'FX160': 57, 'FX250': 58, 'FX350': 59, 'FX500': 60,
            'FWHMSD': 61, 'FWHM160': 64, 'FWHM250': 65, 'FWHM350': 66,
            'FWHM500': 67}


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


def loo_classify(X, key, label):
    prob = np.full(len(label), np.nan)
    for k in np.unique(key):
        out = key == k
        est = mc.LocalLinear(X[~out], label[~out].astype(float))
        prob[out] = est(X[out])
    return prob


def report(label, truth, pred, phi):
    line = ('  %-34s overall %5.1f%%  found %5.1f%%  kept %5.1f%%  |'
            % (label, 100 * np.mean(truth == pred),
               100 * pred[truth].mean(), 100 * (~pred[~truth]).mean()))
    for lo, hi in zip(PHI_BINS[:-1], PHI_BINS[1:]):
        m = (phi >= lo) & (phi < hi)
        line += ' %5.1f' % (100 * np.mean(truth[m] == pred[m])
                            if m.sum() > 20 else np.nan)
    print(line, flush=True)


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
    base = np.column_stack([lg(S['sig'][ok]), lg(S['conc'][ok]),
                            lg(S['cfoot'][ok]), lg(theta)])
    phi = S['cfoot'][ok]
    node = np.asarray(S['node'])[ok]
    key = node[:, 0] * 100000 + node[:, 1] * 1000 + node[:, 2]

    Gr = read_grid(a.grid)
    ckeys, C = read_named(a.catalog, CAT_COLS)
    cat = {k: {n: C[n][i] for n in CAT_COLS} for i, k in enumerate(ckeys)}
    get = lambda n: np.array([cat[k][n] if k in cat else np.nan             # noqa
                              for k in key])
    xi = np.array([Gr[k]['xi_max'] if k in Gr else np.nan for k in key])

    sizes = np.column_stack([
        lg(get('FWHM500') / get('FWHM250')),
        lg(get('FWHM160') / get('FWHM250')),
        lg(get('FWHMSD') / get('FWHM250'))])
    colours = np.column_stack([
        lg(get('FX160') / get('FX250')),
        lg(get('FX250') / get('FX500'))])

    good = (np.all(np.isfinite(base), axis=1) & np.isfinite(xi)
            & np.all(np.isfinite(sizes), axis=1)
            & np.all(np.isfinite(colours), axis=1))
    base, sizes, colours = base[good], sizes[good], colours[good]
    phi, key, xi = phi[good], key[good], xi[good]
    truth = xi > XI_CRIT
    print('samples %d on %d nodes; truly super-critical %.1f%%\n'
          % (len(key), len(np.unique(key)), 100 * truth.mean()), flush=True)
    print('Classification accuracy, node withheld, probability above 0.5.\n'
          'Columns after the bar: accuracy in bins of the measured footprint '
          'factor phi =\n%s.  Real catalogues have median phi between 1.83 and '
          '2.08.\n'
          % ', '.join('%.1f-%.1f' % (l, h)
                      for l, h in zip(PHI_BINS[:-1], PHI_BINS[1:])),
          flush=True)

    for lab, X in (('four observables', base),
                   ('+ per-band size ratios', np.hstack([base, sizes])),
                   ('+ spectral colours', np.hstack([base, colours])),
                   ('+ both', np.hstack([base, sizes, colours])),
                   ('sizes and colours only',
                    np.hstack([sizes, colours]))):
        report(lab, truth, loo_classify(X, key, truth) > 0.5, phi)


if __name__ == '__main__':
    main()
