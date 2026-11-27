#!/usr/bin/env python3
"""
dump_lookup_tables.py -- write out the lookup tables the corrections use, so
that they can be inspected directly.

The corrections are lookups in a table of models.  The table is assembled in
memory at start-up by build_recovery_samples() and never written to disk, which
makes it invisible.  This script writes it out, and also writes, for a chosen
query point, the neighbouring rows that the lookup actually used and the answer
they gave.  Reading those two files together shows exactly what happens when a
real source is corrected: which models it resembled, what their correction
factors were, and how the trend through them was evaluated.

OUTPUTS
    lookup_mass.txt        the full table: one row per model measured through
                           one footprint, with the four observables, the
                           correction factor, the reported and true masses, and
                           the grid indices of the model
    lookup_mass_query.txt  for each query, the LL_NEIGHBOURS nearest rows in the
                           order the lookup weights them, with the distance, the
                           weight, and each row's correction factor, followed by
                           the value the local fit returned

The distance is measured after each observable has been divided by its own
standard deviation over the table, which is what makes the four comparable.
The weight is the tricube function of the distance divided by the distance to
the furthest neighbour used, so the nearest rows count most and the furthest
count nothing.

USAGE
    python3 dump_lookup_tables.py
    python3 dump_lookup_tables.py --query 5.0e21 2.5 1.84 30.0
    python3 dump_lookup_tables.py --max-rows 2000      (a fragment, for reading)
"""
import argparse
import sys

sys.path.insert(0, 'python')
import numpy as np                                     # noqa: E402
import mass_correction_pipeline_6 as mc                # noqa: E402

lg = np.log10


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--catalog', default='cats/bes_model_grid_final3_catalog')
    ap.add_argument('--rectab',
                    default='cats/bes_model_grid_final3_recovery_tables')
    ap.add_argument('--out', default='lookup_mass.txt')
    ap.add_argument('--max-rows', type=int, default=0,
                    help='write only this many rows, evenly spaced, so the '
                         'file can be read by eye')
    ap.add_argument('--query', type=float, nargs=4, action='append',
                    metavar=('SIGMA', 'ZETA', 'PHI', 'RFOOT'),
                    help='a source to look up: background surface density in '
                         'cm^-2, concentration, footprint factor, footprint '
                         'radius in arcsec. May be repeated.')
    a = ap.parse_args()
    mc.VERBOSE = 0

    S = mc.build_recovery_samples(a.catalog, a.rectab)
    ok = np.isfinite(S['rbe_pc']) & (S['rbe_pc'] > 0)
    theta = S['rbe_pc'][ok] * 206264.806 / mc.GRID_DISTANCE_PC
    sig, conc = S['sig'][ok], S['conc'][ok]
    phi, cm = S['cfoot'][ok], S['cm'][ok]
    mrec, mbe = S['mrec'][ok], S['mbe'][ok]
    node = np.asarray(S['node'])[ok]
    X = np.column_stack([lg(sig), lg(conc), lg(phi), lg(theta)])
    good = np.all(np.isfinite(X), axis=1) & np.isfinite(cm) & (cm > 0)
    X, sig, conc, phi, theta = (X[good], sig[good], conc[good], phi[good],
                                theta[good])
    cm, mrec, mbe, node = cm[good], mrec[good], mbe[good], node[good]

    idx = np.arange(len(cm))
    if a.max_rows and a.max_rows < len(cm):
        idx = np.linspace(0, len(cm) - 1, a.max_rows).astype(int)

    with open(a.out, 'w') as fh:
        fh.write('# Lookup table of the mass correction.  One row per model '
                 'measured through one\n# footprint.  A real source is '
                 'corrected by finding the rows whose first four\n# columns '
                 'most resemble its own and reading off C_M.\n#\n')
        fh.write('# %d rows from %d models; %s of them written here\n'
                 % (len(cm), len(np.unique(node[:, 0] * 100000
                                           + node[:, 1] * 1000 + node[:, 2])),
                    'all' if len(idx) == len(cm) else str(len(idx))))
        fh.write('#\n#  1 SIGMA     cm^-2   background surface density\n'
                 '#  2 ZETA       -       concentration, peak over mean in '
                 'the footprint\n'
                 '#  3 PHI        -       footprint factor, footprint over '
                 'source width\n'
                 '#  4 RFOOT      arcsec  footprint radius\n'
                 '#  5 C_M        -       correction factor, true mass over '
                 'measured mass\n'
                 '#  6 M_MEAS     Msun    mass the measurement would report\n'
                 '#  7 M_TRUE     Msun    true Bonnor-Ebert mass of the model\n'
                 '#  8 i  9 j 10 k  -     grid indices of the model\n#\n')
        fh.write('#%12s %11s %9s %9s %10s %12s %12s %5s %4s %4s\n'
                 % ('SIGMA', 'ZETA', 'PHI', 'RFOOT', 'C_M', 'M_MEAS',
                    'M_TRUE', 'i', 'j', 'k'))
        for i in idx:
            fh.write('%13.5e %11.4f %9.4f %9.2f %10.4f %12.5e %12.5e '
                     '%5d %4d %4d\n'
                     % (sig[i], conc[i], phi[i], theta[i], cm[i], mrec[i],
                        mbe[i], node[i, 0], node[i, 1], node[i, 2]))
    print('wrote %s (%d rows)' % (a.out, len(idx)))

    queries = a.query or [[5.0e21, 2.5, 1.84, 30.0], [1.0e22, 4.0, 2.05, 55.0]]
    scale = X.std(axis=0)
    est = mc.LocalLinear(X, lg(cm))
    qout = a.out.replace('.txt', '_query.txt')
    with open(qout, 'w') as fh:
        fh.write('# The rows the lookup actually used, for each query source.\n'
                 '# DIST is measured after each observable is divided by its '
                 'own standard\n# deviation over the table; WEIGHT is the '
                 'tricube weight the fit gave that row.\n')
        for q in queries:
            qv = np.array([lg(q[0]), lg(q[1]), lg(q[2]), lg(q[3])])
            d = np.sqrt((((X - qv) / scale) ** 2).sum(axis=1))
            o = np.argsort(d)[:mc.LL_NEIGHBOURS]
            dmax = d[o[-1]]
            w = (1 - (d[o] / dmax) ** 3) ** 3
            ans = float(10 ** est(qv.reshape(1, -1))[0])
            fh.write('\n# QUERY: SIGMA=%.4e  ZETA=%.4f  PHI=%.4f  '
                     'RFOOT=%.2f arcsec\n' % tuple(q))
            fh.write('# neighbours used: %d; answer returned by the local fit: '
                     'C_M = %.4f\n' % (len(o), ans))
            fh.write('# simple weighted mean of the neighbours, for comparison:'
                     ' C_M = %.4f\n'
                     % float(10 ** (np.sum(w * lg(cm[o])) / np.sum(w))))
            fh.write('#%9s %8s %13s %11s %9s %9s %10s %5s %4s %4s\n'
                     % ('DIST', 'WEIGHT', 'SIGMA', 'ZETA', 'PHI', 'RFOOT',
                        'C_M', 'i', 'j', 'k'))
            for rank, i in enumerate(o):
                fh.write('%10.4f %8.4f %13.5e %11.4f %9.4f %9.2f %10.4f '
                         '%5d %4d %4d\n'
                         % (d[i], w[rank], sig[i], conc[i], phi[i], theta[i],
                            cm[i], node[i, 0], node[i, 1], node[i, 2]))
    print('wrote %s (%d queries)' % (qout, len(queries)))


if __name__ == '__main__':
    main()
