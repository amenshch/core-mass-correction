#!/usr/bin/env python3
"""
chromatic_v2.py -- the stability classification with per-band widths measured
the same way on both sides.

WHAT CHANGED SINCE THE FIRST ATTEMPT
------------------------------------
The first test took the model-side per-band sizes from the 67-column catalogue,
FWHM160bsl and its companions, which are the widths of the whole model on a
noiseless image.  The extraction side supplied getsf's AFWHM, which is the width
of the extracted source after getsf's own background subtraction and within its
own footprint.  These are different measurements, the ratio of the 500 to the
250 micrometre width being 1.118 on one side and 1.279 on the other, and a
classifier trained on one and applied to the other fails.

Here the model side instead uses fwhm_rec from the recovery table: the width of
the truncated, rim-subtracted profile, tabulated per band against the same
truncation radius as the concentration and the footprint factor, which are the
observables already shown to transfer correctly to real extractions.  The
comparison is therefore between quantities constructed in the same way.

Three width ratios are formed, chosen to keep as many samples as possible given
that the width is undefined where a band is unresolved at that truncation:
250 to the surface-density band, 350 to 250, and 500 to 250.  Their getsf
counterparts in the injection catalogues are AFWHM04/AFWHM03, AFWHM05/AFWHM04
and AFWHM06/AFWHM04, the band assignment being 01 for 70, 02 for 160, 03 the
surface density, 04 for 250, 05 for 350 and 06 for 500 micrometres.

USAGE
    python3 chromatic_v2.py               grid test and injection test
    python3 chromatic_v2.py --test grid   grid test only
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
PHI_BINS = [1.0, 1.4, 1.7, 1.9, 2.1, 2.5, 4.0]
RUNS = ['260811_cmf2p00_n600_s2030', '260812_cmf2p00_n600_s2031',
        '260813_cmf2p00_n600_s2032', '260814_cmf2p00_n600_s2033']
VARIANTS = ['', '_flat', '_iso', '_flat_iso']


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


def parse_truth(path):
    out = []
    for ln in open(path):
        if not ln.strip() or ln.lstrip()[0] in '#!':
            continue
        f = ln.split()
        t = f[1]
        out.append((float(f[2]), float(f[3]),
                    (int(t[1:3]), int(t[4:6]), int(t[7:9]))))
    return out


def report(label, truth, pred, phi=None):
    line = ('  %-38s accuracy %5.1f%%  found %5.1f%%  kept %5.1f%%'
            % (label, 100 * np.mean(truth == pred),
               100 * pred[truth].mean() if truth.any() else np.nan,
               100 * (~pred[~truth]).mean() if (~truth).any() else np.nan))
    if phi is not None:
        line += '  |'
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
    ap.add_argument('--grid2', default='cats/bes_grid_final2.txt')
    ap.add_argument('--injdir', default='Aquila-injections')
    ap.add_argument('--test', default='both')
    ap.add_argument('--with500', action='store_true',
                    help='include the 500 to 250 micrometre width ratio, whose '
                         'model and getsf values disagree because the 500 '
                         'micrometre width is beam-determined')
    ap.add_argument('--variants', nargs='+', default=VARIANTS,
                    help="injection variants to use; '' is as placed")
    a = ap.parse_args()
    mc.VERBOSE = 0

    S = mc.build_recovery_samples(a.catalog, a.rectab)
    ok = np.isfinite(S['rbe_pc']) & (S['rbe_pc'] > 0)
    theta = S['rbe_pc'][ok] * 206264.806 / mc.GRID_DISTANCE_PC
    base = np.column_stack([lg(S['sig'][ok]), lg(S['conc'][ok]),
                            lg(S['cfoot'][ok]), lg(theta)])
    ratios = [('250/SD', S['w250'][ok] / S['w_sd'][ok]),
              ('350/250', S['w350'][ok] / S['w250'][ok])]
    if a.with500:
        ratios.append(('500/250', S['w500'][ok] / S['w250'][ok]))
    wid = np.column_stack([lg(v) for _, v in ratios])
    phi = S['cfoot'][ok]
    node = np.asarray(S['node'])[ok]
    key = node[:, 0] * 100000 + node[:, 1] * 1000 + node[:, 2]
    Gr = read_grid(a.grid)
    xi = np.array([Gr[(n[0], n[1], n[2])]['xi_max']
                   if (n[0], n[1], n[2]) in Gr else np.nan for n in node])
    g = (np.all(np.isfinite(base), axis=1) & np.all(np.isfinite(wid), axis=1)
         & np.isfinite(xi))
    base, wid, phi, key, xi = base[g], wid[g], phi[g], key[g], xi[g]
    label = xi > XI_CRIT
    print('training samples with all width ratios finite: %d on %d nodes, '
          '%.1f%% super-critical\n'
          % (len(key), len(np.unique(key)), 100 * label.mean()), flush=True)
    print('median model width ratios: %s'
          % ', '.join('%s %.3f' % (n, m) for (n, _), m
                      in zip(ratios, np.median(10 ** wid, axis=0))),
          flush=True)

    if a.test in ('grid', 'both'):
        print('\nGrid test, leave-one-node-out.  Columns after the bar: '
              'accuracy in bins of the\nmeasured footprint factor phi = %s.\n'
              % ', '.join('%.1f-%.1f' % (l, h)
                          for l, h in zip(PHI_BINS[:-1], PHI_BINS[1:])),
              flush=True)
        for lab, X in (('four observables', base),
                       ('four plus profile width ratios',
                        np.hstack([base, wid]))):
            p = np.full(len(key), np.nan)
            for k in np.unique(key):
                o = key == k
                p[o] = mc.LocalLinear(X[~o], label[~o].astype(float))(X[o])
            report(lab, label, p > 0.5, phi)

    if a.test in ('inject', 'both'):
        grid2 = read_grid(a.grid2)
        rows = []
        for r in RUNS:
            for v in a.variants:
                d = os.path.join(a.injdir, r + v)
                cats = [os.path.join(d, 'Aquila.s.sources.ok.cat'),
                        os.path.join(d, 'Aquila.s.sources.ok.add.cat')]
                if not os.path.exists(cats[0]):
                    continue
                table, _ = mc.correct_cloud(cats, DIST)
                t = GetsfTable(cats[0])
                col = lambda n: np.asarray(t.col(n),                        # noqa
                                           float).reshape(t.nrows, -1)[:, 0]
                idx = {int(n): i for i, n in enumerate(col('NO'))}
                xs, ys = col('XCO_P'), col('YCO_P')
                aSD, a250 = col('AFWHM03'), col('AFWHM04')
                a350, a500 = col('AFWHM05'), col('AFWHM06')
                tr = parse_truth(glob.glob(os.path.join(d, '*_truth.txt'))[0])
                tx = np.array([q[0] for q in tr])
                ty = np.array([q[1] for q in tr])
                for i, no in enumerate(table['NO']):
                    j = idx.get(int(no))
                    if j is None:
                        continue
                    dd = np.hypot(tx - xs[j], ty - ys[j])
                    m = int(np.argmin(dd))
                    if dd[m] >= MATCH_PIX or tr[m][2] not in grid2:
                        continue
                    rows.append((table['sigma'][i], table['conc'][i],
                                 table['rfoot'][i], table['foot_as'][i],
                                 a250[j] / aSD[j], a350[j] / a250[j],
                                 a500[j] / a250[j],
                                 grid2[tr[m][2]]['xi_max'],
                                 VARIANTS.index(v)))
        D = np.array(rows)
        qb = np.column_stack([lg(D[:, 0]), lg(D[:, 1]), lg(D[:, 2]),
                              lg(D[:, 3])])
        cols = [4, 5] + ([6] if a.with500 else [])
        qw = np.column_stack([lg(D[:, c]) for c in cols])
        ts = D[:, 7] > XI_CRIT
        var = D[:, 8].astype(int)
        v = np.all(np.isfinite(qb), axis=1) & np.all(np.isfinite(qw), axis=1)
        qb, qw, ts, var = qb[v], qw[v], ts[v], var[v]
        print('\nInjection test: %d matched sources, %.1f%% truly '
              'super-critical, so the commoner\nclass alone would score %.1f%%.'
              % (len(ts), 100 * ts.mean(),
                 100 * max(ts.mean(), 1 - ts.mean())), flush=True)
        print('  median getsf width ratios: %s'
              % ', '.join('%s %.3f' % (n, m) for (n, _), m
                          in zip(ratios, np.median(10 ** qw, axis=0))),
              flush=True)
        print('  NOTE: each injected core appears once per variant, so pooled '
              'entries are not\n  independent; the per-variant blocks below '
              'are.\n', flush=True)
        preds = {}
        for lab, Xtr, Xq in (('four observables', base, qb),
                             ('four plus width ratios',
                              np.hstack([base, wid]), np.hstack([qb, qw]))):
            preds[lab] = mc.LocalLinear(Xtr, label.astype(float))(Xq) > 0.5
        for vi, vname in [(None, 'pooled over the four variants')] + \
                [(VARIANTS.index(x), x or 'as placed') for x in a.variants]:
            m = np.ones(len(ts), bool) if vi is None else (var == vi)
            if m.sum() < 20:
                continue
            print('  %s: %d sources, %d truly super-critical'
                  % (vname, int(m.sum()), int(ts[m].sum())), flush=True)
            for lab in preds:
                sens = (100 * preds[lab][m][ts[m]].mean() if ts[m].any()
                        else np.nan)
                spec = (100 * (~preds[lab][m][~ts[m]]).mean()
                        if (~ts[m]).any() else np.nan)
                print('     %-32s found %5.1f%%  kept %5.1f%%  balanced %5.1f%%'
                      % (lab, sens, spec, 0.5 * (sens + spec)), flush=True)


if __name__ == '__main__':
    main()
