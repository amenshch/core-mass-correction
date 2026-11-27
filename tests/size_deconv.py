"""Can the model grid replace Gaussian deconvolution for half-maximum sizes?

The observer measures a width at half maximum that is broadened by the beam and
removes the beam by subtracting it in quadrature.  The quantity sought is the
intrinsic width at half maximum of the source.  In the grid that quantity is
known exactly: it is the axis FWHM_pc of the model, converted to arcseconds at
the distance at which the grid was rendered.

Estimators compared, all from quantities an observer measures:
    half-maximum width as measured, w_sd, undeconvolved;
    the Gaussian deconvolution sqrt(w_sd^2 - beam^2), standard practice;
    a grid inversion trained on (Sigma, zeta, phi, footprint radius, w_sd).

The grid inversion is validated by withholding one node at a time, then applied
to sources injected into real Aquila maps and re-extracted by getsf, whose
intrinsic widths are known from the grid node each was drawn from.  The
extraction supplies AFWHM on the surface-density band in place of w_sd, so the
injection test also measures whether the inversion transfers, which is where the
per-band width ratios failed earlier.
"""
import glob, os, sys
sys.path.insert(0, 'python'); sys.path.insert(0, '.')
import numpy as np
import mass_correction_pipeline_6 as mc
from getsf_columns import GetsfTable
lg = np.log10
BEAM, DIST, MATCH = 13.5, 260.0, 10.0
RUNS = ['260811_cmf2p00_n600_s2030', '260812_cmf2p00_n600_s2031',
        '260813_cmf2p00_n600_s2032', '260814_cmf2p00_n600_s2033']
VAR = {'': 'as placed', '_flat': 'flattened', '_iso': 'isolated',
       '_flat_iso': 'flattened and isolated'}


def read_grid(p):
    rows = [l.split() for l in open(p) if l.strip()]
    c = {n: i for i, n in enumerate(rows[0])}
    out = {}
    for r in rows[1:]:
        d = {}
        for n, i in c.items():
            try: d[n] = float(r[i])
            except ValueError: d[n] = np.nan
        out[(int(r[c['i']]), int(r[c['j']]), int(r[c['k']]))] = d
    return out


def q(label, ratio):
    g = np.isfinite(ratio) & (ratio > 0)
    if g.sum() < 5:
        print('     %-40s too few' % label); return
    lr = lg(ratio[g])
    print('     %-40s n=%5d  median %6.3f  scatter %5.2f  within 1.5x %3.0f%%'
          '  typical error %6.3f dex'
          % (label, int(g.sum()), np.median(ratio[g]), 10 ** lr.std(),
             100 * np.mean((ratio[g] > 1 / 1.5) & (ratio[g] < 1.5)),
             np.median(np.abs(lr))), flush=True)


mc.VERBOSE = 0
S = mc.build_recovery_samples('cats/bes_model_grid_final3_catalog',
                              'cats/bes_model_grid_final3_recovery_tables')
ok = np.isfinite(S['rbe_pc']) & (S['rbe_pc'] > 0)
w = S['w_sd'][ok]; phi = S['cfoot'][ok]; r_as = phi * w / 2.0
node = np.asarray(S['node'])[ok]
key = node[:, 0] * 100000 + node[:, 1] * 1000 + node[:, 2]
G3 = read_grid('cats/bes_grid_final3.txt')
SC = 206264.806 / mc.GRID_DISTANCE_PC
fw_true = np.array([G3[tuple(n)]['FWHM_pc'] * SC if tuple(n) in G3 else np.nan
                    for n in node])
X = np.column_stack([lg(S['sig'][ok]), lg(S['conc'][ok]), lg(phi), lg(r_as),
                     lg(w)])
g = np.all(np.isfinite(X), axis=1) & np.isfinite(fw_true) & (fw_true > 0)
X, key, fw_true, w, r_as = X[g], key[g], fw_true[g], w[g], r_as[g]
print('grid: %d samples on %d nodes; true intrinsic FWHM median %.1f arcsec\n'
      % (len(key), len(np.unique(key)), np.median(fw_true)), flush=True)

pred = np.full(len(key), np.nan)
for k in np.unique(key):
    o = key == k
    pred[o] = mc.LocalLinear(X[~o], lg(fw_true[~o]))(X[o])
pred = 10 ** pred
dec = np.sqrt(np.clip(w ** 2 - BEAM ** 2, 0, None))
print('Grid test, ratio of each estimate of the intrinsic half-maximum width '
      'to its true value:')
q('measured width, undeconvolved', w / fw_true)
q('Gaussian deconvolution', dec / fw_true)
q('grid inversion, node withheld', pred / fw_true)

est = mc.LocalLinear(X, lg(fw_true))
G2 = read_grid('cats/bes_grid_final2.txt')
print('\nInjection test, the same estimators applied to getsf-extracted '
      'sources.', flush=True)
for suf, name in VAR.items():
    rows = []
    for r in RUNS:
        d = os.path.join('Aquila-injections', r + suf)
        cats = [os.path.join(d, 'Aquila.s.sources.ok.cat'),
                os.path.join(d, 'Aquila.s.sources.ok.add.cat')]
        if not os.path.exists(cats[0]):
            continue
        tab, _ = mc.correct_cloud(cats, DIST)
        t = GetsfTable(cats[0])
        cc = lambda n: np.asarray(t.col(n), float).reshape(t.nrows, -1)[:, 0]
        idx = {int(n): i for i, n in enumerate(cc('NO'))}
        xs, ys, aSD = cc('XCO_P'), cc('YCO_P'), cc('AFWHM03')
        tr = []
        for ln in open(glob.glob(os.path.join(d, '*_truth.txt'))[0]):
            if not ln.strip() or ln.lstrip()[0] in '#!':
                continue
            f = ln.split(); tg = f[1]
            tr.append((float(f[2]), float(f[3]),
                       (int(tg[1:3]), int(tg[4:6]), int(tg[7:9]))))
        tx = np.array([a[0] for a in tr]); ty = np.array([a[1] for a in tr])
        for i, no in enumerate(tab['NO']):
            j = idx.get(int(no))
            if j is None:
                continue
            dd = np.hypot(tx - xs[j], ty - ys[j]); k = int(np.argmin(dd))
            if dd[k] >= MATCH or tr[k][2] not in G2:
                continue
            rows.append((tab['sigma'][i], tab['conc'][i], tab['rfoot'][i],
                         tab['foot_as'][i], aSD[j],
                         G2[tr[k][2]]['FWHM_pc'] * SC))
    D = np.array(rows)
    if len(D) == 0:
        continue
    qq = np.column_stack([lg(D[:, 0]), lg(D[:, 1]), lg(D[:, 2]), lg(D[:, 3]),
                          lg(D[:, 4])])
    v = np.all(np.isfinite(qq), axis=1) & np.isfinite(D[:, 5]) & (D[:, 5] > 0)
    qq, ws, ft = qq[v], D[v, 4], D[v, 5]
    print('\n  %s: %d sources, true intrinsic FWHM median %.1f arcsec'
          % (name, len(ft), np.median(ft)), flush=True)
    q('measured width, undeconvolved', ws / ft)
    q('Gaussian deconvolution', np.sqrt(np.clip(ws ** 2 - BEAM ** 2, 0,
                                                None)) / ft)
    q('grid inversion', 10 ** est(qq) / ft)
