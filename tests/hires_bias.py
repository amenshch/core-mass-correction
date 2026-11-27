#!/usr/bin/env python3
"""
hires_bias.py -- how much does the corrector move when its observables are
measured on a hires-derived surface-density map instead of the true one?

The correctors are trained on quantities read from the model's true column
density, and applied to quantities read from a map that hires reconstructs by
fitting a modified blackbody in every pixel of the 160-500 micrometre images.
The two are not the same measurement.  This script quantifies the difference
from paired azimuthally averaged radial profiles of the same model stamp, one
from the true surface density and one from the hires reconstruction, both at
13.5 arcsec resolution on 3 arcsec pixels.

Everything is computed with the pipeline's own definitions, since those are what
the corrector consumes: at a truncation radius r the rim value is the profile
value at r, the profile is truncated and rim-subtracted there, fwhm_rec is the
width at half maximum of that truncated profile, the footprint factor is
phi = 2 r / fwhm_rec, and the concentration zeta is the rim-subtracted peak
divided by the area-weighted rim-subtracted mean inside r.  An asymptotic cloud
level must NOT be used in place of the rim value: it removes the cancellation
that rim subtraction provides and greatly overstates the effect.

The two maps are compared at matched measured phi, not at matched truncation
radius, because phi is what an observer measures and what the corrector is
queried with.

The last points of each profile are discarded: azimuthal averaging over the
corners of a square stamp drops the mean there.

USAGE
    python3 hires_bias.py [--phi 1.84 2.10] [--dir hires]
"""
import argparse
import os
import sys

sys.path.insert(0, 'python')
sys.path.insert(0, '.')
import numpy as np                                     # noqa: E402
import mass_correction_pipeline_6 as mc                # noqa: E402

lg = np.log10
NODES = ['i02j15k12', 'i02j23k10', 'i02j25k09', 'i05j15k11', 'i05j23k10',
         'i05j25k09', 'i08j17k11', 'i08j23k11', 'i08j25k09', 'i05j16k01',
         'i05j24k12', 'i05j27k16', 'i04j16k13', 'i08j23k12', 'i05j30k14',
         'i09j22k13']


def read_profile(path):
    r, v = [], []
    for line in open(path, encoding='utf-8', errors='ignore'):
        if line.startswith('#') or not line.strip():
            continue
        f = line.split()
        r.append(float(f[1])); v.append(float(f[3]))
    r, v = np.array(r), np.array(v)
    good = np.isfinite(v) & (v > 1e-30)
    return r[good], v[good]


def trim_corners(r, v, frac=0.85):
    """Drop the outermost points, where azimuthal averaging over the corners of
    a square stamp depresses the mean.  Inspection of the profiles shows the
    drop setting in at about 0.89 of the maximum radius in every node, so 0.85
    is a safe cut; 0.90 is not."""
    keep = r <= frac * r.max()
    return r[keep], v[keep]


def measure(x, y, rt):
    rim = float(np.interp(rt, x, y))
    m = x <= rt
    xx, s = x[m], y[m] - rim
    if len(xx) < 4 or s[0] <= 0:
        return None
    half = s[0] / 2.0
    fw = np.nan
    for i in range(len(xx) - 1):
        if s[i] >= half > s[i + 1]:
            fw = 2 * (xx[i] + (xx[i + 1] - xx[i]) * (s[i] - half)
                      / (s[i] - s[i + 1]))
            break
    if not np.isfinite(fw) or fw <= 0:
        return None
    mean = np.trapezoid(s * 2 * np.pi * xx, xx) / (np.pi * rt ** 2)
    if mean <= 0:
        return None
    return dict(rt=rt, fwhm=fw, phi=2 * rt / fw, zeta=s[0] / mean)


def curve(x, y):
    """Measurements over a ladder of truncation radii, as a function of phi."""
    out = []
    for rt in np.arange(x[1], 0.8 * x.max(), (x[1] - x[0])):
        m = measure(x, y, rt)
        if m:
            out.append(m)
    if not out:
        return None
    phi = np.array([o['phi'] for o in out])
    o = np.argsort(phi)
    return dict(phi=phi[o], zeta=np.array([q['zeta'] for q in out])[o],
                fwhm=np.array([q['fwhm'] for q in out])[o],
                rt=np.array([q['rt'] for q in out])[o])


def at(c, phi):
    if c is None or phi < c['phi'].min() or phi > c['phi'].max():
        return None
    return {k: float(np.interp(phi, c['phi'], c[k]))
            for k in ('zeta', 'fwhm', 'rt')}


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--dir', default='hires')
    ap.add_argument('--phi', type=float, nargs='+', default=[1.84, 2.10])
    ap.add_argument('--catalog', default='cats/bes_model_grid_final3_catalog')
    ap.add_argument('--rectab',
                    default='cats/bes_model_grid_final3_recovery_tables')
    a = ap.parse_args()
    mc.VERBOSE = 0

    S = mc.build_recovery_samples(a.catalog, a.rectab)
    ok = np.isfinite(S['rbe_pc']) & (S['rbe_pc'] > 0)
    th = S['rbe_pc'][ok] * 206264.806 / mc.GRID_DISTANCE_PC
    X = np.column_stack([lg(S['sig'][ok]), lg(S['conc'][ok]),
                         lg(S['cfoot'][ok]), lg(th)])
    cm = S['cm'][ok]
    g = np.all(np.isfinite(X), axis=1) & np.isfinite(cm) & (cm > 0)
    est = mc.LocalLinear(X[g], lg(cm[g]))

    print('Background level of each map, and the observables compared at '
          'matched measured phi.\nAll ratios are hires divided by true.\n',
          flush=True)
    hdr = ('%-11s %9s |' % ('node', 'bg ratio'))
    for p in a.phi:
        hdr += '  phi=%.2f: %6s %6s %7s' % (p, 'zeta', 'fwhm', 'dC_M')
    print(hdr, flush=True)
    rows = []
    for n in NODES:
        pn = os.path.join(a.dir, n + '_nc.dat')
        ph = os.path.join(a.dir, n + '_hi.dat')
        if not (os.path.exists(pn) and os.path.exists(ph)):
            print('%-11s  profiles missing' % n)
            continue
        rn0, vn0 = read_profile(pn)
        rh0, vh0 = read_profile(ph)
        # Keep only radii present and positive in both maps, so the two
        # profiles are compared on one common grid.
        common = np.intersect1d(rn0, rh0)
        vn0 = vn0[np.isin(rn0, common)]
        vh0 = vh0[np.isin(rh0, common)]
        rn, vn = trim_corners(common, vn0)
        rh, vh = trim_corners(common, vh0)
        # Background from the flat outer part of each profile, on the same
        # radial window for both maps so the ratio is not distorted.
        far = rn >= 0.6 * rn.max()
        bgn = float(np.median(vn[far]))
        bgh = float(np.median(vh[far]))
        cn, ch = curve(rn, vn), curve(rh, vh)
        line = '%-11s %9.4f |' % (n, bgh / bgn)
        rec = dict(node=n, bg=bgh / bgn)
        for p in a.phi:
            an, ah = at(cn, p), at(ch, p)
            if an is None or ah is None:
                line += '  %-6s %6s %6s %7s' % ('', '-', '-', '-')
                continue
            q0 = np.array([[lg(bgn), lg(an['zeta']), lg(p), lg(an['rt'])]])
            q1 = np.array([[lg(bgh), lg(ah['zeta']), lg(p), lg(ah['rt'])]])
            c0 = float(10 ** est(q0)[0]); c1 = float(10 ** est(q1)[0])
            line += '           %6.3f %6.3f %6.1f%%' % (
                ah['zeta'] / an['zeta'], ah['fwhm'] / an['fwhm'],
                100 * (c1 / c0 - 1))
            rec['d%.2f' % p] = 100 * (c1 / c0 - 1)
        rows.append(rec)
        print(line, flush=True)

    print('\nsummary over the nodes measured:')
    print('  background ratio: median %.4f, range %.4f to %.4f'
          % (np.median([r['bg'] for r in rows]),
             min(r['bg'] for r in rows), max(r['bg'] for r in rows)))
    for p in a.phi:
        v = [r['d%.2f' % p] for r in rows if 'd%.2f' % p in r]
        if v:
            print('  change in the correction factor at phi=%.2f: median '
                  '%+.1f%%, range %+.1f%% to %+.1f%% over %d nodes'
                  % (p, np.median(v), min(v), max(v), len(v)))


if __name__ == '__main__':
    main()
