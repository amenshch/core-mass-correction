#!/usr/bin/env python3
"""
invert_grid.py -- inversion of the four measured observables onto the physical
parameters of the embedded Bonnor-Ebert models, and its validation.

PURPOSE
-------
The mass correction uses the model grid to answer one question, by how much a
reported mass must be multiplied.  The grid contains much more: every model is
generated from a known physical state, and knows its own dimensionless
truncation radius xi_max, its central-to-edge density contrast, its true mass
and its stability.  This script asks whether that state can be recovered from
the same four observables the corrector uses, namely

    Sigma  the background surface density on the surface-density band,
    zeta   the peak-to-mean concentration,
    phi    the footprint-to-size ratio,
    theta  the angle subtended by the truncation radius, in arcseconds,

which are the only ones already validated end to end against real extractions.

The script has two modes.  In validation mode it withholds one grid node at a
time, infers the physical parameters of every sample of that node from the
remaining nodes, and compares them with the truth.  Nothing should be inferred
about real cores until this shows the inversion to be determinate: if the four
observables do not fix a parameter within the grid itself, no statement about
that parameter for a real core is meaningful.  In application mode it inverts
the real cores of the corrected per-cloud catalogs.

ESTIMATOR
---------
The same local linear regression the corrector uses, LocalLinear from the
pipeline, with the same neighbour count, applied to the base-ten logarithm of
the target parameter.  Using the corrector's own estimator means that any
weakness found here is a property of the observables and the grid, not of a
newly introduced fitting method.

USAGE
    python3 invert_grid.py --validate
    python3 invert_grid.py --apply --clouds Aquila OrionA
    python3 invert_grid.py --validate --targets xi_max contr_rho M_BE

Inputs expected in the current directory:
    python/mass_correction_pipeline_6.py, python/getsf_columns.py
    cats/bes_model_grid_final3_catalog
    cats/bes_model_grid_final3_recovery_tables
    cats/bes_grid_final3.txt
    out_f3_snr1.5/<cloud>.corrected.cat            (application mode only)

QUANTITIES REPORTED IN VALIDATION MODE, per target parameter:
    n                number of grid samples scored
    typical error    median of the absolute base-ten logarithm of the ratio of
                     the inferred to the true value, in decimal exponents
    scatter          standard deviation of that logarithm, as a factor
    within 1.5x      percentage of samples inferred within a factor of 1.5
    class accuracy   for xi_max only: percentage of samples whose stability
                     class, sub-critical or super-critical with respect to
                     xi_max = 6.451, is correctly recovered
"""
import argparse
import os
import sys

sys.path.insert(0, 'python')
import numpy as np                                     # noqa: E402
import mass_correction_pipeline_6 as mc                # noqa: E402

lg = np.log10
XI_CRIT = 6.451          # dimensionless truncation radius of the critical sphere
GRIDFILE = 'cats/bes_grid_final3.txt'
CATALOG = 'cats/bes_model_grid_final3_catalog'
RECTAB = 'cats/bes_model_grid_final3_recovery_tables'
CLOUDS = ['Scorpius', 'Ophiuchus', 'Aquila', 'OrionA', 'California',
          'CygnusX', 'W3W4W5']
TARGETS_DEFAULT = ['xi_max', 'contr_rho', 'M_BE', 'rho_c']


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
                d[n] = np.nan      # non-numeric columns, e.g. the stability label
        out[key] = d
    return out


def build(catalog, rectab, gridfile, targets):
    """Return the four-observable design matrix, the node key, and the truth
    value of each target for every training sample."""
    S = mc.build_recovery_samples(catalog, rectab)
    ok = np.isfinite(S['rbe_pc']) & (S['rbe_pc'] > 0)
    theta = S['rbe_pc'][ok] * 206264.806 / mc.GRID_DISTANCE_PC
    X = np.column_stack([lg(S['sig'][ok]), lg(S['conc'][ok]),
                         lg(S['cfoot'][ok]), lg(theta)])
    node = np.asarray(S['node'])[ok]
    key = node[:, 0] * 100000 + node[:, 1] * 1000 + node[:, 2]
    G = read_grid(gridfile)
    truth = {}
    for t in targets:
        truth[t] = np.array([G[k][t] if k in G else np.nan for k in key])
    good = np.all(np.isfinite(X), axis=1)
    for t in targets:
        good &= np.isfinite(truth[t]) & (truth[t] > 0)
    return X[good], key[good], {t: truth[t][good] for t in targets}


def validate(X, key, truth, targets):
    uniq = np.unique(key)
    print('samples = %d, nodes = %d' % (len(key), len(uniq)))
    pred = {t: np.full(len(key), np.nan) for t in targets}
    for kk in uniq:
        out = key == kk
        keep = ~out
        for t in targets:
            est = mc.LocalLinear(X[keep], lg(truth[t][keep]))
            pred[t][out] = est(X[out])
    print('\n%-12s %8s %14s %9s %12s %15s'
          % ('parameter', 'n', 'typical error', 'scatter', 'within 1.5x',
             'class accuracy'))
    for t in targets:
        r = 10 ** pred[t] / truth[t]
        g = np.isfinite(r) & (r > 0)
        lr = lg(r[g])
        acc = ''
        if t == 'xi_max':
            c_true = truth[t][g] > XI_CRIT
            c_pred = 10 ** pred[t][g] > XI_CRIT
            acc = '%13.1f%%' % (100 * np.mean(c_true == c_pred))
            base = max(c_true.mean(), 1 - c_true.mean())
            print('  stability classification with respect to xi_max = %.3f:'
                  % XI_CRIT)
            print('    truly super-critical %.1f%% of samples, so a classifier '
                  'that always answered' % (100 * c_true.mean()))
            print('    with the commoner class would score %.1f%%' % (100 * base))
            for lab, sel in (('truly super-critical', c_true),
                             ('truly sub-critical', ~c_true)):
                if sel.sum():
                    print('    %-22s n=%6d, called super-critical %.1f%% of '
                          'the time' % (lab, sel.sum(),
                                        100 * c_pred[sel].mean()))
        print('%-12s %8d %13.3f %9.2f %10.0f%% %15s'
              % (t, g.sum(), np.median(np.abs(lr)), 10 ** lr.std(),
                 100 * np.mean((r[g] > 1 / 1.5) & (r[g] < 1.5)), acc))
    return pred


def apply_to_clouds(X, truth, targets, clouds, outdir):
    est = {t: mc.LocalLinear(X, lg(truth[t])) for t in targets}
    print('\n%-11s %7s %10s %12s %14s %14s'
          % ('cloud', 'n', 'in hull', 'median xi_max', 'super-critical',
             'median contrast'))
    for name in clouds:
        p = os.path.join(outdir, '%s.corrected.cat' % name)
        if not os.path.exists(p):
            print('%-11s  no corrected catalog at %s' % (name, p))
            continue
        rows = [l.split() for l in open(p) if l.strip() and l[0] not in '!#']
        a = np.array([[float(r[3]), float(r[4]), float(r[5]), float(r[26]),
                       float(r[32])] for r in rows])
        q = np.column_stack([lg(a[:, 0]), lg(a[:, 1]), lg(a[:, 2]),
                             lg(a[:, 3])])
        inh = a[:, 4] > 0.5
        g = np.all(np.isfinite(q), axis=1) & inh
        xi = 10 ** est['xi_max'](q[g])
        cr = 10 ** est['contr_rho'](q[g])
        print('%-11s %7d %9.1f%% %12.2f %13.1f%% %14.1f'
              % (name, int(g.sum()), 100 * inh.mean(), np.median(xi),
                 100 * np.mean(xi > XI_CRIT), np.median(cr)))


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--validate', action='store_true')
    ap.add_argument('--apply', action='store_true')
    ap.add_argument('--catalog', default=CATALOG)
    ap.add_argument('--rectab', default=RECTAB)
    ap.add_argument('--grid', default=GRIDFILE)
    ap.add_argument('--outdir', default='out_f3_snr1.5')
    ap.add_argument('--clouds', nargs='+', default=CLOUDS)
    ap.add_argument('--targets', nargs='+', default=TARGETS_DEFAULT)
    a = ap.parse_args()
    if not (a.validate or a.apply):
        a.validate = True
    mc.VERBOSE = 0
    X, key, truth = build(a.catalog, a.rectab, a.grid, a.targets)
    if a.validate:
        validate(X, key, truth, a.targets)
    if a.apply:
        apply_to_clouds(X, truth, a.targets, a.clouds, a.outdir)


if __name__ == '__main__':
    main()
