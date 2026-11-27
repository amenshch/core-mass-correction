#!/usr/bin/env python3
"""
loo_neighbour_scan.py -- leave-one-out error of the mass correction as a
function of the number of training samples used per local-linear prediction.

MOTIVATION
----------
The estimator used by the pipeline, mass_correction_pipeline_6.LocalLinear, is
a tricube-weighted local linear regression that fits a plane to the
LL_NEIGHBOURS nearest training samples of the query point, distances being
measured after each axis has been divided by its own standard deviation.  The
constant LL_NEIGHBOURS = 320 was chosen on the final2 model grid, where it was
4.3 percent of the 7497 training samples.  On the extended final3 grid the same
constant is 1.1 percent of the 30 022 samples, so the volume of the
four-dimensional neighborhood over which the plane is fitted has contracted by
about a factor of three, and the linear extent of that neighborhood by about a
factor of 1.4 in each coordinate.  The final2 and final3 leave-one-out numbers
are therefore not produced by an identically smoothed estimator.  This script
measures the dependence directly, so that the neighbor count can be chosen on
evidence rather than inherited.

METHOD
------
Identical to loo_4axis_full_v3.py in every respect except that the neighbor
count is varied.  One node is withheld at a time; both correctors are built on
the remaining nodes; the withheld node is scored at the single training sample
whose footprint factor is closest to REP_RFOOT, the value typical of the grid.
The three-observable corrector uses the base-ten logarithms of the surface
density, the concentration, and the footprint factor; the four-observable
corrector adds the base-ten logarithm of the angle subtended by the truncation
radius, in arcseconds.

USAGE
    Put in the current directory:
        mass_correction_pipeline_6.py
        getsf_columns.py
        bes_model_grid_final3_catalog
        bes_model_grid_final3_recovery_tables
        bes_grid_final3.txt
    then:
        python3 loo_neighbour_scan.py
        python3 loo_neighbour_scan.py --k 160 320 640 --out scan.csv

QUANTITIES REPORTED, for each neighbor count and each corrector:
    med ratio      median over the nodes of the ratio of the corrected mass to
                   the true Bonnor-Ebert mass M_BE of the node
    scatter        standard deviation of the base-ten logarithm of that ratio,
                   expressed as a multiplicative factor
    within 2x      percentage of nodes with the ratio between 0.5 and 2
    typ err        median of the absolute base-ten logarithm of the ratio, in
                   decimal exponents, over all scored nodes
    typ err, hole  the same, restricted to nodes carrying the hole flag, that
                   is, nodes whose density at the outer truncation boundary
                   falls below the density of the embedding cloud
    bounded        number of nodes at which the local plane extrapolated beyond
                   the range of its own neighbors and the prediction had to be
                   clipped to that range
"""
import argparse
import csv
import os
import sys
import time

sys.path.insert(0, '.')
import numpy as np                                     # noqa: E402
import mass_correction_pipeline_6 as mc                # noqa: E402

lg = np.log10
CATALOG = 'bes_model_grid_final3_catalog'
RECTAB = 'bes_model_grid_final3_recovery_tables'
GRIDFILE = 'bes_grid_final3.txt'
REP_RFOOT = 1.84
K_DEFAULT = (80, 160, 320, 640, 1280, 2560)


def read_hole(path):
    rows = [l.split() for l in open(path) if l.strip()]
    col = {n: i for i, n in enumerate(rows[0])}
    return {int(r[col['i']]) * 100000 + int(r[col['j']]) * 1000
            + int(r[col['k']]): int(r[col['hole']]) for r in rows[1:]}


def summarize(ratio):
    lr = lg(ratio)
    return (np.median(ratio), 10 ** lr.std(),
            100 * np.mean((ratio > 0.5) & (ratio < 2)), np.median(np.abs(lr)))


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--catalog', default=CATALOG)
    ap.add_argument('--rectab', default=RECTAB)
    ap.add_argument('--grid', default=GRIDFILE)
    ap.add_argument('--k', type=int, nargs='+', default=list(K_DEFAULT))
    ap.add_argument('--out', default='loo_neighbour_scan.csv')
    a = ap.parse_args()

    for f in (a.catalog, a.rectab, a.grid, 'mass_correction_pipeline_6.py'):
        if not os.path.exists(f):
            sys.exit('ERROR: required input file not found: %s' % f)

    hole = read_hole(a.grid)
    S = mc.build_recovery_samples(a.catalog, a.rectab)
    ok = np.isfinite(S['rbe_pc']) & (S['rbe_pc'] > 0)
    node = np.asarray(S['node'])[ok]
    theta = S['rbe_pc'][ok] * 206264.806 / mc.GRID_DISTANCE_PC
    sig, conc, cfoot = S['sig'][ok], S['conc'][ok], S['cfoot'][ok]
    cm, mrec, mbe = S['cm'][ok], S['mrec'][ok], S['mbe'][ok]

    X3 = np.column_stack([lg(sig), lg(conc), lg(cfoot)])
    X4 = np.column_stack([lg(sig), lg(conc), lg(cfoot), lg(theta)])
    y = lg(cm)
    key = node[:, 0] * 100000 + node[:, 1] * 1000 + node[:, 2]
    uniq = np.unique(key)
    print('samples = %d, nodes = %d, neighbor counts to scan: %s'
          % (len(key), len(uniq), a.k))

    rep = [np.where(key == kk)[0][
        np.argmin(np.abs(cfoot[np.where(key == kk)[0]] - REP_RFOOT))]
        for kk in uniq]
    is_hole = np.array([hole.get(int(kk), 0) == 1 for kk in uniq])
    print('nodes carrying the hole flag: %d of %d'
          % (int(is_hole.sum()), len(uniq)))

    rows = []
    for K in a.k:
        t0 = time.time()
        p3 = np.full(len(uniq), np.nan)
        p4 = np.full(len(uniq), np.nan)
        b3 = np.zeros(len(uniq), bool)
        b4 = np.zeros(len(uniq), bool)
        for i, kk in enumerate(uniq):
            out = key == kk
            keep = ~out
            j = rep[i]
            e3 = mc.LocalLinear(X3[keep], y[keep], k=K)
            e4 = mc.LocalLinear(X4[keep], y[keep], k=K)
            v3, v4 = e3(X3[j:j + 1]), e4(X4[j:j + 1])
            p3[i], p4[i] = float(v3[0]), float(v4[0])
            b3[i], b4[i] = bool(e3._flag[0]), bool(e4._flag[0])
        good = np.isfinite(p3) & np.isfinite(p4)
        m_be = mbe[rep]
        m_rec = mrec[rep]
        r3 = (m_rec * 10 ** p3)[good] / m_be[good]
        r4 = (m_rec * 10 ** p4)[good] / m_be[good]
        h = is_hole[good]
        for lab, r, bb in (('three observables', r3, b3),
                           ('four observables', r4, b4)):
            med, sc, w2, te = summarize(r)
            teh = np.median(np.abs(lg(r[h]))) if h.sum() else np.nan
            rows.append(dict(k=K, corrector=lab, n=int(good.sum()),
                             med_ratio=med, scatter=sc, within2x=w2,
                             typ_err=te, typ_err_hole=teh,
                             bounded=int(bb[good].sum())))
        print('  K = %5d done in %5.0f s' % (K, time.time() - t0))

    print('\n%-6s %-20s %5s %9s %8s %9s %9s %11s %8s'
          % ('K', 'corrector', 'n', 'med ratio', 'scatter', 'within2x',
             'typ err', 'typ err hole', 'bounded'))
    for r in rows:
        print('%-6d %-20s %5d %9.3f %8.2f %8.0f%% %9.4f %11.4f %8d'
              % (r['k'], r['corrector'], r['n'], r['med_ratio'], r['scatter'],
                 r['within2x'], r['typ_err'], r['typ_err_hole'], r['bounded']))
    with open(a.out, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        for r in rows:
            w.writerow(r)
    print('\nWritten %s' % a.out)


if __name__ == '__main__':
    main()
