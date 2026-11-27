#!/usr/bin/env python3
"""
loo_compare_grids.py -- side-by-side summary of two leave-one-out runs of
loo_4axis_full_v3.py, one on the final2 model grid and one on the extended
final3 model grid.

WHY A SEPARATE SCRIPT
---------------------
The final3 grid extends to embedding background surface densities a factor of
eight lower, and to Bonnor-Ebert masses a factor of five lower, than final2.
A leave-one-out error computed over the whole final3 volume is therefore
measured over a larger and intrinsically harder region of parameter space, and
is not directly comparable with the final2 number.  This script splits the
final3 nodes into those lying inside the parameter volume spanned by final2 and
those lying outside it, so that the like-for-like comparison can be read off
separately from the behavior in the newly added region.  The nodes are also
split by the hole flag, a hole node being one whose density at the outer
truncation boundary falls below the density of the embedding cloud.

INPUT FILES, all in the current directory:
    loo4d_full_results_final2_rerun.csv   output of loo_4axis_full_v3.py on final2
    loo4d_full_results_final3.csv         output of loo_4axis_full_v3.py on final3
    bes_grid_final2.txt                   final2 grid parameter file
    bes_grid_final3.txt                   final3 grid parameter file

USAGE
    python3 loo_compare_grids.py
    python3 loo_compare_grids.py --csv2 FILE --csv3 FILE --grid2 FILE --grid3 FILE

QUANTITIES REPORTED, in each block and for each of the three mass estimates
(the uncorrected reported mass, the mass corrected on three observables, and
the mass corrected on four observables):
    n              number of grid nodes contributing to the block
    median         median over the nodes of the ratio of the mass obtained to
                   the true Bonnor-Ebert mass M_BE of the node
    scatter        standard deviation of the base-ten logarithm of that ratio,
                   expressed as a multiplicative factor
    within 2x      percentage of nodes whose ratio lies between 0.5 and 2
    med|log10|     median of the absolute value of the base-ten logarithm of
                   the ratio, that is, the typical per-node error
The final line of each block gives the paired per-node comparison: the number
of nodes at which the four-observable corrector is closer to M_BE than the
three-observable corrector, and the number at which it is farther.
"""
import argparse
import csv
import sys

import numpy as np

QUANT = (('sd', 'SD_emb'), ('m', 'M_BE'), ('fw', 'FWHM_pc'))


def read_grid(path):
    """Return {node_key: parameter dict} from a bes_grid_*.txt file.

    The node key is i * 10**5 + j * 10**3 + k, the same key written by
    loo_4axis_full_v3.py.
    """
    rows = [l.split() for l in open(path) if l.strip()]
    col = {n: i for i, n in enumerate(rows[0])}
    out = {}
    for r in rows[1:]:
        key = (int(r[col['i']]) * 100000 + int(r[col['j']]) * 1000
               + int(r[col['k']]))
        d = {short: float(r[col[name]]) for short, name in QUANT}
        d['hole'] = int(r[col['hole']])
        out[key] = d
    return out


def read_loo(path):
    """Return the array of (node_key, M_BE, M_reported, M_3axis, M_4axis)."""
    rec = []
    for r in csv.DictReader(open(path)):
        if not r['M_corrected_3axis_Msun'] or not r['M_corrected_4axis_Msun']:
            continue
        rec.append((int(r['node_key']), float(r['M_BE_Msun']),
                    float(r['M_reported_Msun']),
                    float(r['M_corrected_3axis_Msun']),
                    float(r['M_corrected_4axis_Msun'])))
    return rec


def report(label, rec):
    if not rec:
        print('  %s: no nodes' % label)
        return
    a = np.array([x[1:] for x in rec])
    print('  %s   (n = %d nodes)' % (label, len(rec)))
    print('     %-34s %8s %9s %10s %11s'
          % ('', 'median', 'scatter', 'within 2x', 'med|log10|'))
    for name, j in (('reported mass, uncorrected', 1),
                    ('corrected, three observables', 2),
                    ('corrected, four observables', 3)):
        r = a[:, j] / a[:, 0]
        lr = np.log10(r)
        print('     %-34s %8.3f %9.2f %9.0f%% %11.3f'
              % (name, np.median(r), 10 ** lr.std(),
                 100 * np.mean((r > 0.5) & (r < 2)), np.median(np.abs(lr))))
    d3 = np.abs(np.log10(a[:, 2] / a[:, 0]))
    d4 = np.abs(np.log10(a[:, 3] / a[:, 0]))
    print('     four observables closer to M_BE at %d nodes, farther at %d'
          % (int((d4 < d3).sum()), int((d4 > d3).sum())))


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--csv2', default='loo4d_full_results_final2_rerun.csv')
    ap.add_argument('--csv3', default='loo4d_full_results_final3.csv')
    ap.add_argument('--grid2', default='bes_grid_final2.txt')
    ap.add_argument('--grid3', default='bes_grid_final3.txt')
    a = ap.parse_args()

    for f in (a.csv2, a.csv3, a.grid2, a.grid3):
        try:
            open(f).close()
        except OSError:
            sys.exit('ERROR: required input file not found: %s' % f)

    g2, g3 = read_grid(a.grid2), read_grid(a.grid3)
    r2, r3 = read_loo(a.csv2), read_loo(a.csv3)

    lim = {q: (min(v[q] for v in g2.values()), max(v[q] for v in g2.values()))
           for q, _ in QUANT}
    print('Parameter volume spanned by the final2 grid, used below to define '
          'the restricted comparison:')
    for q, name in QUANT:
        print('   %-9s %10.4g  to %10.4g' % (name, lim[q][0], lim[q][1]))

    def inside(key):
        p = g3[key]
        return all(lim[q][0] * 0.999 <= p[q] <= lim[q][1] * 1.001
                   for q, _ in QUANT)

    print('\nA. final2 grid, all nodes')
    report('final2', r2)
    print('\nB. final3 grid, all nodes')
    report('final3, complete', r3)
    print('\nC. final3 grid, nodes inside the final2 parameter volume '
          '(the like-for-like comparison with block A)')
    report('final3, restricted', [x for x in r3 if inside(x[0])])
    print('\nD. final3 grid, nodes outside the final2 parameter volume '
          '(the newly added region)')
    report('final3, new region', [x for x in r3 if not inside(x[0])])
    print('\nE. final3 grid, nodes with the hole flag not set')
    report('final3, hole = 0', [x for x in r3 if g3[x[0]]['hole'] == 0])
    print('\nF. final3 grid, nodes with the hole flag set')
    report('final3, hole = 1', [x for x in r3 if g3[x[0]]['hole'] == 1])


if __name__ == '__main__':
    main()
