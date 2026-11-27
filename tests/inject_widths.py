#!/usr/bin/env python3
"""
inject_widths.py -- do the per-band profile width ratios improve the mass
correction and the temperature correction on real extractions?

BACKGROUND
----------
On the model grid, adding the ratios of the per-band widths of the truncated,
rim-subtracted profile to the four validated observables reduces the
leave-one-node-out error of the mass factor C_M = M_BE / M_ex from 0.028 to
0.015 dex, and of the temperature factor C_T = Td_avg / T_SED from 0.003 to
0.001 dex.  The same addition, however, degrades the sub- against
super-critical classification when applied to getsf-extracted sources, because
a systematic offset of a few per cent between the model-side and
extraction-side widths is enough to move sources across a decision boundary.

A regression on a smoothly varying factor need not be as fragile as a decision
boundary, so the question is open and is settled here directly: both correctors
are trained on the model grid and applied to sources injected into real Aquila
maps and re-extracted by getsf, whose true masses and true mass-averaged dust
temperatures are known.

WHAT IS COMPARED
----------------
For the mass, the ratio of the obtained mass to the injected model's true
Bonnor-Ebert mass, for three treatments: the reported mass uncorrected, the
mass corrected on the four observables, and the mass corrected on the four
observables plus the two width ratios.

For the temperature, the ratio of the obtained temperature to the injected
model's mass-averaged dust temperature Td_avg, for two treatments: the
temperature from the fit to the spectral energy distribution as reported, and
that temperature multiplied by the predicted C_T.

Each injection variant is reported separately, since a given core appears once
in each variant and pooled entries are therefore not independent.

USAGE
    python3 inject_widths.py
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
DIST = 260.0
MATCH_PIX = 10.0
RUNS = ['260811_cmf2p00_n600_s2030', '260812_cmf2p00_n600_s2031',
        '260813_cmf2p00_n600_s2032', '260814_cmf2p00_n600_s2033']
VARIANTS = ['', '_flat', '_iso', '_flat_iso']
VNAME = {'': 'as placed', '_flat': 'flattened', '_iso': 'isolated',
         '_flat_iso': 'flattened and isolated'}
# 1-based columns of the 67-column model catalogue
COL_TDAVG, COL_TSED4, COL_TSED3 = 12, 17, 23


def read_cat_cols(path, cols):
    keys, out = [], {n: [] for n in cols}
    for line in open(path):
        if line.startswith('#') or not line.strip():
            continue
        f = line.split()
        keys.append((int(f[1]), int(f[2]), int(f[3])))
        for n, c in cols.items():
            try:
                out[n].append(float(f[c - 1]))
            except (ValueError, IndexError):
                out[n].append(np.nan)
    return {k: {n: out[n][i] for n in cols} for i, k in enumerate(keys)}


def parse_truth(path):
    out = []
    for ln in open(path):
        if not ln.strip() or ln.lstrip()[0] in '#!':
            continue
        f = ln.split()
        t = f[1]
        out.append((float(f[2]), float(f[3]),
                    (int(t[1:3]), int(t[4:6]), int(t[7:9])), float(f[8])))
    return out


def stats(label, ratio):
    g = np.isfinite(ratio) & (ratio > 0)
    if g.sum() < 5:
        print('     %-34s too few sources' % label)
        return
    lr = lg(ratio[g])
    print('     %-34s n=%3d  median %6.3f  scatter %5.2f  within 2x %3.0f%%  '
          'typical error %6.3f dex'
          % (label, int(g.sum()), np.median(ratio[g]), 10 ** lr.std(),
             100 * np.mean((ratio[g] > 0.5) & (ratio[g] < 2)),
             np.median(np.abs(lr))), flush=True)


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--catalog', default='cats/bes_model_grid_final3_catalog')
    ap.add_argument('--rectab',
                    default='cats/bes_model_grid_final3_recovery_tables')
    ap.add_argument('--cat2', default='cats/bes_model_grid_final2_catalog')
    ap.add_argument('--injdir', default='Aquila-injections')
    a = ap.parse_args()
    mc.VERBOSE = 0

    # ---- training on the extended grid -------------------------------------
    S = mc.build_recovery_samples(a.catalog, a.rectab)
    ok = np.isfinite(S['rbe_pc']) & (S['rbe_pc'] > 0)
    theta = S['rbe_pc'][ok] * 206264.806 / mc.GRID_DISTANCE_PC
    base = np.column_stack([lg(S['sig'][ok]), lg(S['conc'][ok]),
                            lg(S['cfoot'][ok]), lg(theta)])
    wid = np.column_stack([lg(S['w250'][ok] / S['w_sd'][ok]),
                           lg(S['w350'][ok] / S['w250'][ok])])
    cm = S['cm'][ok]
    node = np.asarray(S['node'])[ok]
    cat3 = read_cat_cols(a.catalog, {'Td': COL_TDAVG, 'T4': COL_TSED4,
                                     'T3': COL_TSED3})
    td = np.array([cat3[tuple(n)]['Td'] if tuple(n) in cat3 else np.nan
                   for n in node])
    t4 = np.array([cat3[tuple(n)]['T4'] if tuple(n) in cat3 else np.nan
                   for n in node])
    t3 = np.array([cat3[tuple(n)]['T3'] if tuple(n) in cat3 else np.nan
                   for n in node])
    tsed = np.where(np.isfinite(t4) & (t4 > 0), t4, t3)
    ct = td / tsed

    fin = (np.all(np.isfinite(base), axis=1) & np.isfinite(cm) & (cm > 0)
           & np.isfinite(ct) & (ct > 0))
    finw = fin & np.all(np.isfinite(wid), axis=1)
    print('training: %d samples for the four-observable correctors, %d for the '
          'width-augmented ones\n' % (int(fin.sum()), int(finw.sum())),
          flush=True)
    est = {
        'M4': mc.LocalLinear(base[fin], lg(cm[fin])),
        'M4w': mc.LocalLinear(np.hstack([base, wid])[finw], lg(cm[finw])),
        'T4': mc.LocalLinear(base[fin], lg(ct[fin])),
        'T4w': mc.LocalLinear(np.hstack([base, wid])[finw], lg(ct[finw])),
    }

    # ---- injected sources ---------------------------------------------------
    cat2 = read_cat_cols(a.cat2, {'Td': COL_TDAVG})
    per = {v: [] for v in VARIANTS}
    for r in RUNS:
        for v in VARIANTS:
            d = os.path.join(a.injdir, r + v)
            cats = [os.path.join(d, 'Aquila.s.sources.ok.cat'),
                    os.path.join(d, 'Aquila.s.sources.ok.add.cat')]
            if not os.path.exists(cats[0]):
                continue
            table, _ = mc.correct_cloud(cats, DIST)
            t = GetsfTable(cats[0])
            col = lambda n: np.asarray(t.col(n),                            # noqa
                                       float).reshape(t.nrows, -1)[:, 0]
            idx = {int(n): i for i, n in enumerate(col('NO'))}
            xs, ys = col('XCO_P'), col('YCO_P')
            aSD, a250, a350 = col('AFWHM03'), col('AFWHM04'), col('AFWHM05')
            tr = parse_truth(glob.glob(os.path.join(d, '*_truth.txt'))[0])
            tx = np.array([q[0] for q in tr])
            ty = np.array([q[1] for q in tr])
            for i, no in enumerate(table['NO']):
                j = idx.get(int(no))
                if j is None:
                    continue
                dd = np.hypot(tx - xs[j], ty - ys[j])
                m = int(np.argmin(dd))
                if dd[m] >= MATCH_PIX:
                    continue
                ijk = tr[m][2]
                per[v].append((table['sigma'][i], table['conc'][i],
                               table['rfoot'][i], table['foot_as'][i],
                               a250[j] / aSD[j], a350[j] / a250[j],
                               table['M_init'][i], tr[m][3],
                               table['T_init'][i],
                               cat2[ijk]['Td'] if ijk in cat2 else np.nan))

    print('\nMass: ratio of the obtained mass to the injected model\'s true '
          'Bonnor-Ebert mass.\nTemperature: ratio of the obtained temperature '
          'to the model\'s mass-averaged dust\ntemperature Td_avg.  Variants '
          'are independent of one another only across cores,\nnot within them.',
          flush=True)
    for v in VARIANTS:
        D = np.array(per[v])
        if len(D) == 0:
            continue
        qb = np.column_stack([lg(D[:, 0]), lg(D[:, 1]), lg(D[:, 2]),
                              lg(D[:, 3])])
        qw = np.column_stack([lg(D[:, 4]), lg(D[:, 5])])
        mex, mbe, tin, tdt = D[:, 6], D[:, 7], D[:, 8], D[:, 9]
        g = (np.all(np.isfinite(qb), axis=1) & np.all(np.isfinite(qw), axis=1)
             & np.isfinite(mex) & (mex > 0) & np.isfinite(mbe) & (mbe > 0))
        qb, qw, mex, mbe, tin, tdt = (qb[g], qw[g], mex[g], mbe[g], tin[g],
                                      tdt[g])
        qbw = np.hstack([qb, qw])
        print('\n  %s: %d matched sources' % (VNAME[v], len(mex)), flush=True)
        stats('mass, uncorrected', mex / mbe)
        stats('mass, four observables', mex * 10 ** est['M4'](qb) / mbe)
        stats('mass, four plus width ratios',
              mex * 10 ** est['M4w'](qbw) / mbe)
        stats('temperature, as fitted', tin / tdt)
        stats('temperature, four observables',
              tin * 10 ** est['T4'](qb) / tdt)
        stats('temperature, four plus width ratios',
              tin * 10 ** est['T4w'](qbw) / tdt)


if __name__ == '__main__':
    main()
