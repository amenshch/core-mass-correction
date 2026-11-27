#!/usr/bin/env python3
"""
inject_chromatic.py -- does the gain from the per-band size ratios survive real
extraction, and do those ratios help the other inverted quantities?

TEST A: THE STABILITY CLASSIFICATION ON INJECTED SOURCES
--------------------------------------------------------
On the model grid, adding the ratios of the sizes measured in different bands
to the four validated observables raises the accuracy of the sub-critical
against super-critical classification from 82--86 per cent to about 96.5 per
cent at the footprint factors real catalogues reach.  Those model sizes are
measured on noiseless images over the whole source.  Here the same classifier,
trained on the grid, is applied to sources injected into real Aquila maps and
re-extracted by getsf, whose per-band sizes carry the noise, the
background-subtraction differences between bands, and the beam limitation at
500 micrometres that the model sizes do not.  The truth is the injected model's
own dimensionless truncation radius.  The comparison of interest is between the
four observables alone and the four plus the size ratios, on the same sources.

Note that the injected models are drawn from the earlier grid and are absent
from the extended grid on which the classifier is trained, so this is an
out-of-training test in the strict sense.

Band assignment in the injection catalogues, read from their headers: 01 is 70,
02 is 160, 03 the surface density, 04 is 250, 05 is 350 and 06 is 500
micrometres.  The size ratios are formed to match those used on the model side:
500 to 250, 160 to 250, and surface density to 250.

TEST B: THE OTHER INVERTED QUANTITIES
-------------------------------------
The same additional observables are given to the leave-one-node-out inversion
of the mass factor C_M, the temperature factor C_T and the truncation radius
xi_max, to see whether the improvement is confined to the classification.

USAGE
    python3 inject_chromatic.py                 both tests
    python3 inject_chromatic.py --test A        injected sources only
"""
import argparse
import glob
import os
import sys

sys.path.insert(0, 'python')
sys.path.insert(0, '.')
import numpy as np                                     # noqa: E402
import mass_correction_pipeline_6 as mc                # noqa: E402
from getsf_columns import GetsfTable                   # noqa: E402

lg = np.log10
XI_CRIT = 6.451
DIST = 260.0
MATCH_PIX = 10.0
CAT_COLS = {'FWHMSD': 61, 'FWHM160': 64, 'FWHM250': 65, 'FWHM500': 67}
RUNS = ['260811_cmf2p00_n600_s2030', '260812_cmf2p00_n600_s2031',
        '260813_cmf2p00_n600_s2032', '260814_cmf2p00_n600_s2033']
VARIANTS = ['', '_flat', '_iso', '_flat_iso']


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
        d = {}
        for n, i in col.items():
            try:
                d[n] = float(r[i])
            except ValueError:
                d[n] = np.nan
        out[(int(r[col['i']]), int(r[col['j']]), int(r[col['k']]))] = d
    return out


def build_training(catalog, rectab, gridfile):
    S = mc.build_recovery_samples(catalog, rectab)
    ok = np.isfinite(S['rbe_pc']) & (S['rbe_pc'] > 0)
    theta = S['rbe_pc'][ok] * 206264.806 / mc.GRID_DISTANCE_PC
    base = np.column_stack([lg(S['sig'][ok]), lg(S['conc'][ok]),
                            lg(S['cfoot'][ok]), lg(theta)])
    node = np.asarray(S['node'])[ok]
    key = node[:, 0] * 100000 + node[:, 1] * 1000 + node[:, 2]
    Gr = read_grid(gridfile)
    xi = np.array([Gr[(n[0], n[1], n[2])]['xi_max']
                   if (n[0], n[1], n[2]) in Gr else np.nan for n in node])
    ckeys, C = read_named(catalog, CAT_COLS)
    cat = {k: {n: C[n][i] for n in CAT_COLS} for i, k in enumerate(ckeys)}
    g = lambda n: np.array([cat[k][n] if k in cat else np.nan for k in key])  # noqa
    sizes = np.column_stack([lg(g('FWHM500') / g('FWHM250')),
                             lg(g('FWHM160') / g('FWHM250')),
                             lg(g('FWHMSD') / g('FWHM250'))])
    good = (np.all(np.isfinite(base), axis=1) & np.isfinite(xi)
            & np.all(np.isfinite(sizes), axis=1))
    return (base[good], sizes[good], key[good], xi[good],
            S, ok, np.asarray(S['cm'])[ok][good])


def parse_truth(path):
    out = []
    for ln in open(path):
        if not ln.strip() or ln.lstrip()[0] in '#!':
            continue
        f = ln.split()
        tag = f[1]
        ijk = (int(tag[1:3]), int(tag[4:6]), int(tag[7:9]))
        out.append(dict(x=float(f[2]), y=float(f[3]), ijk=ijk))
    return out


def load_injection(run_dir, grid2):
    """Return per matched source: four observables, size ratios, true xi_max."""
    cats = [os.path.join(run_dir, 'Aquila.s.sources.ok.cat'),
            os.path.join(run_dir, 'Aquila.s.sources.ok.add.cat')]
    if not os.path.exists(cats[0]):
        return None
    table, _ = mc.correct_cloud(cats, DIST)
    t = GetsfTable(cats[0])
    col = lambda n: np.asarray(t.col(n), float).reshape(t.nrows, -1)[:, 0]  # noqa
    nos = col('NO')
    idx = {int(n): i for i, n in enumerate(nos)}
    xs, ys = col('XCO_P'), col('YCO_P')
    a160, aSD, a250, a500 = (col('AFWHM02'), col('AFWHM03'), col('AFWHM04'),
                             col('AFWHM06'))
    truth = parse_truth(glob.glob(os.path.join(run_dir, '*_truth.txt'))[0])
    tx = np.array([r['x'] for r in truth])
    ty = np.array([r['y'] for r in truth])
    rows = []
    for i, no in enumerate(table['NO']):
        j = idx.get(int(no))
        if j is None:
            continue
        d = np.hypot(tx - xs[j], ty - ys[j])
        m = int(np.argmin(d))
        if d[m] >= MATCH_PIX:
            continue
        ijk = truth[m]['ijk']
        if ijk not in grid2:
            continue
        rows.append((table['sigma'][i], table['conc'][i], table['rfoot'][i],
                     table['foot_as'][i], a500[j] / a250[j],
                     a160[j] / a250[j], aSD[j] / a250[j],
                     grid2[ijk]['xi_max']))
    return np.array(rows) if rows else None


def loo_class(X, key, label):
    p = np.full(len(label), np.nan)
    for k in np.unique(key):
        o = key == k
        p[o] = mc.LocalLinear(X[~o], label[~o].astype(float))(X[o])
    return p


def loo_reg(X, key, y):
    p = np.full(len(y), np.nan)
    for k in np.unique(key):
        o = key == k
        p[o] = mc.LocalLinear(X[~o], lg(y[~o]))(X[o])
    return 10 ** p


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--catalog', default='cats/bes_model_grid_final3_catalog')
    ap.add_argument('--rectab',
                    default='cats/bes_model_grid_final3_recovery_tables')
    ap.add_argument('--grid', default='cats/bes_grid_final3.txt')
    ap.add_argument('--grid2', default='cats/bes_grid_final2.txt')
    ap.add_argument('--injdir', default='Aquila-injections')
    ap.add_argument('--test', default='AB')
    a = ap.parse_args()
    mc.VERBOSE = -1

    base, sizes, key, xi, S, ok, cm = build_training(a.catalog, a.rectab,
                                                     a.grid)
    label = xi > XI_CRIT
    print('training: %d samples on %d nodes, %.1f%% super-critical\n'
          % (len(key), len(np.unique(key)), 100 * label.mean()), flush=True)

    if 'A' in a.test:
        grid2 = read_grid(a.grid2)
        rows = []
        for r in RUNS:
            for v in VARIANTS:
                out = load_injection(os.path.join(a.injdir, r + v), grid2)
                if out is not None:
                    rows.append(out)
        D = np.vstack(rows)
        q_base = np.column_stack([lg(D[:, 0]), lg(D[:, 1]), lg(D[:, 2]),
                                  lg(D[:, 3])])
        q_size = np.column_stack([lg(D[:, 4]), lg(D[:, 5]), lg(D[:, 6])])
        tsuper = D[:, 7] > XI_CRIT
        g = (np.all(np.isfinite(q_base), axis=1)
             & np.all(np.isfinite(q_size), axis=1))
        q_base, q_size, tsuper = q_base[g], q_size[g], tsuper[g]
        print('Test A: injected sources re-extracted by getsf, pooled over the '
              'four seeds and\nthe four variants.  %d matched sources with all '
              'observables finite, of which\n%.1f%% are truly super-critical, '
              'so always answering with the commoner class\nwould score %.1f%%.'
              % (len(tsuper), 100 * tsuper.mean(),
                 100 * max(tsuper.mean(), 1 - tsuper.mean())), flush=True)
        for lab, Xtr, Xq in (('four observables', base, q_base),
                             ('four observables plus size ratios',
                              np.hstack([base, sizes]),
                              np.hstack([q_base, q_size]))):
            est = mc.LocalLinear(Xtr, label.astype(float))
            pred = est(Xq) > 0.5
            print('  %-36s accuracy %5.1f%%  super-critical found %5.1f%%  '
                  'sub-critical kept %5.1f%%'
                  % (lab, 100 * np.mean(pred == tsuper),
                     100 * pred[tsuper].mean(),
                     100 * (~pred[~tsuper]).mean()), flush=True)
        print('  median ratio of the getsf 500 to 250 micrometre sizes: %.3f; '
              'on the model side\n  the same ratio has median %.3f'
              % (np.median(10 ** q_size[:, 0]),
                 np.median(10 ** sizes[:, 0])), flush=True)

    if 'B' in a.test:
        print('\nTest B: leave-one-node-out inversion of other quantities, '
              'with and without the\nper-band size ratios.  The error is the '
              'median absolute base-ten logarithm of\nthe ratio of the '
              'inferred to the true value, in decimal exponents.\n', flush=True)
        ckeys, C = read_named(a.catalog, {'Td_avg': 12, 'T_SED4bsl': 17,
                                          'T_SED3bsl': 23})
        cat = {k: {n: C[n][i] for n in C} for i, k in enumerate(ckeys)}
        gg = lambda n: np.array([cat[k][n] if k in cat else np.nan          # noqa
                                 for k in key])
        ts = gg('T_SED4bsl')
        ts = np.where(np.isfinite(ts) & (ts > 0), ts, gg('T_SED3bsl'))
        targets = [('mass factor C_M', cm), ('truncation radius xi_max', xi),
                   ('temperature factor C_T', gg('Td_avg') / ts)]
        for tlab, y in targets:
            line = '  %-26s' % tlab
            for slab, X in (('four', base),
                            ('plus sizes', np.hstack([base, sizes]))):
                m = np.isfinite(y) & (y > 0)
                pr = loo_reg(X[m], key[m], y[m])
                v = np.isfinite(pr) & (pr > 0)
                line += '  %s %6.3f dex' % (slab,
                                            np.median(np.abs(lg(pr[v]
                                                                / y[m][v]))))
            print(line, flush=True)


if __name__ == '__main__':
    main()
