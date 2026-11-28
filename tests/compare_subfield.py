#!/usr/bin/env python3
"""
compare_subfield.py -- compare two getsf extractions of the same sub-field that
differ only in the 13.5 arcsec surface density band: the original hires image
(uncorrected) and the same image corrected by hicorrect (hicorrected).

Usage (from any folder):
    python3 tests/compare_subfield.py [--dir data/subfield/Aquila] [--dist 260]
                                      [--match 6] [--out FILE]

Inputs, in <dir>/uncorrected and <dir>/hicorrected:
    the getsf catalog merged with its .add.cat (name containing "=" once), the
    same merged with the SED catalog (name containing "=" twice), and the 13.5
    arcsec surface density image (*.fits).
Only sources inside the area covered by BOTH surface density images are
compared, since the two extractions need not cover the same area.

Reported:
    1. detections: sources matched by position (within --match arcsec), and
       sources found in one extraction only;
    2. for matched sources, the ratios corrected/uncorrected of the background
       surface density (PEAK^BGF), the peak surface density above it (PEAK^SRC),
       the concentration zeta, the footprint factor phi, and the angular
       footprint radius theta, all measured in the surface density band;
    3. the mass correction (pipeline version 7, four observables) of the
       sources with QUALITY=ok SED fits: the reported mass M_INIT, the factor
       C_M, and the corrected mass M_CORR, compared source by source.
"""
import argparse
import glob
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                '..', 'method'))
import mass_correction_pipeline_7 as mc                # noqa: E402
import paths                                           # noqa: E402
from astropy.io import fits                            # noqa: E402
from astropy.wcs import WCS                            # noqa: E402

mc.VERBOSE = 0


def pick(d, ncat):
    """File in d whose name contains '=' exactly ncat times."""
    f = [p for p in glob.glob(os.path.join(d, '*=*')) if os.path.basename(p).count('=') == ncat]
    if len(f) != 1:
        sys.exit('expected one catalog with %d "=" in %s, found %s' % (ncat, d, f))
    return f[0]


def coverage(d):
    """Function returning True for (ra, dec) inside the finite, positive area of the
    surface density image in d."""
    f = glob.glob(os.path.join(d, '*.fits'))
    if len(f) != 1:
        sys.exit('expected one FITS image in %s, found %s' % (d, f))
    h = fits.open(f[0])[0]
    img = np.asarray(h.data, float)
    good = np.isfinite(img) & (img > 0)
    w = WCS(h.header).celestial

    def inside(ra, dec):
        x, y = w.world_to_pixel_values(ra, dec)
        i, j = np.round(y).astype(int), np.round(x).astype(int)
        ok = (i >= 0) & (j >= 0) & (i < img.shape[0]) & (j < img.shape[1])
        res = np.zeros(len(ra), bool)
        res[ok] = good[i[ok], j[ok]]
        return res
    return inside


def match(ra1, de1, ra2, de2, rmax):
    """Nearest-neighbor match within rmax arcsec; returns index pairs."""
    c = np.cos(np.radians(np.median(de1)))
    p1, p2 = [], []
    for i in range(len(ra1)):
        d = 3600 * np.hypot((ra2 - ra1[i]) * c, de2 - de1[i])
        j = int(np.argmin(d))
        if d[j] <= rmax:
            p1.append(i); p2.append(j)
    p1, p2 = np.array(p1, int), np.array(p2, int)
    # keep one-to-one pairs only
    _, first = np.unique(p2, return_index=True)
    keep = np.zeros(len(p2), bool); keep[first] = True
    return p1[keep], p2[keep]


def stats(name, r, out):
    r = np.asarray(r, float); r = r[np.isfinite(r) & (r > 0)]
    if len(r) == 0:
        out.append('  %-38s  no values' % name); return
    p10, p50, p90 = np.percentile(r, [10, 50, 90])
    out.append('  %-38s n %4d   median %6.3f   10th %6.3f   90th %6.3f'
               % (name, len(r), p50, p10, p90))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--dir', default=paths.SUBFIELD_DIR)
    ap.add_argument('--dist', type=float, default=260.0)
    ap.add_argument('--match', type=float, default=6.0, help='matching radius, arcsec')
    ap.add_argument('--out', default=None)
    a = ap.parse_args()
    U, C = os.path.join(a.dir, 'uncorrected'), os.path.join(a.dir, 'hicorrected')
    mc.load_getsf._dist = a.dist
    inU, inC = coverage(U), coverage(C)
    out = ['Comparison of two extractions of %s (distance %g pc)' % (a.dir, a.dist),
           'uncorrected: %s' % os.path.basename(pick(U, 1)),
           'hicorrected: %s' % os.path.basename(pick(C, 1)), '']

    # ---------------------------------------------------------- detections
    gU, gC = mc.load_getsf(pick(U, 1)), mc.load_getsf(pick(C, 1))
    comU = inU(gU['ra'], gU['dec']) & inC(gU['ra'], gU['dec'])
    comC = inU(gC['ra'], gC['dec']) & inC(gC['ra'], gC['dec'])
    iU, iC = np.where(comU)[0], np.where(comC)[0]
    mU, mC = match(gU['ra'][iU], gU['dec'][iU], gC['ra'][iC], gC['dec'][iC], a.match)
    pu, pc = iU[mU], iC[mC]
    out += ['1. Detections inside the common area (all sources of the merged catalogs)',
            '  uncorrected %d, hicorrected %d, matched within %g arcsec: %d'
            % (len(iU), len(iC), a.match, len(pu)),
            '  only in uncorrected: %d, only in hicorrected: %d'
            % (len(iU) - len(pu), len(iC) - len(pc)), '']
    # unmatched sources: their peak above background and background, to see who they are
    for lab, g, idx, matched in (('uncorrected', gU, iU, pu), ('hicorrected', gC, iC, pc)):
        un = np.setdiff1d(idx, matched)
        if len(un):
            out.append('  sources only in %s: median PEAK^SRC %.2e, median PEAK^BGF %.2e cm^-2'
                       % (lab, np.nanmedian(g['peak'][un]), np.nanmedian(g['Nbg'][un])))
    out.append('')

    # ------------------------------------------------------- observables
    out.append('2. Ratio hicorrected/uncorrected for matched sources (surface density band)')
    stats('background PEAK^BGF (Sigma)', gC['Nbg'][pc] / gU['Nbg'][pu], out)
    stats('peak above background PEAK^SRC', gC['peak'][pc] / gU['peak'][pu], out)
    stats('contrast (1 + PEAK^SRC/PEAK^BGF)',
          (1 + gC['peak'][pc] / gC['Nbg'][pc]) / (1 + gU['peak'][pu] / gU['Nbg'][pu]), out)
    stats('concentration zeta', gC['conc'][pc] / gU['conc'][pu], out)
    stats('footprint factor phi', gC['rfoot'][pc] / gU['rfoot'][pu], out)
    stats('footprint radius theta', gC['foot_radius_as'][pc] / gU['foot_radius_as'][pu], out)
    out.append('')

    # ---------------------------------------------------- mass correction
    tU, _ = mc.correct_cloud(pick(U, 2), a.dist)
    tC, _ = mc.correct_cloud(pick(C, 2), a.dist)
    tU = {k: np.asarray(v) for k, v in tU.items()}
    tC = {k: np.asarray(v) for k, v in tC.items()}
    kU = inU(tU['ra'], tU['dec']) & inC(tU['ra'], tU['dec'])
    kC = inU(tC['ra'], tC['dec']) & inC(tC['ra'], tC['dec'])
    jU, jC = np.where(kU)[0], np.where(kC)[0]
    nU, nC = match(tU['ra'][jU], tU['dec'][jU], tC['ra'][jC], tC['dec'][jC], a.match)
    qu, qc = jU[nU], jC[nC]
    out += ['3. Mass correction (QUALITY=ok SED fits, common area)',
            '  corrected sources: uncorrected %d, hicorrected %d, matched %d'
            % (len(jU), len(jC), len(qu))]
    for lab, t, j in (('uncorrected', tU, jU), ('hicorrected', tC, jC)):
        out.append('  %s: median C_M %.3f, median C_M_3D %.3f, IN_DOMAIN %d of %d'
                   % (lab, np.nanmedian(t['C_M'][j]), np.nanmedian(t['C_M_3D'][j]),
                      int(np.sum(t['IN_DOMAIN'][j])), len(j)))
    out.append('  ratios hicorrected/uncorrected over matched sources:')
    stats('reported mass M_INIT', tC['M_init'][qc] / tU['M_init'][qu], out)
    stats('SED temperature T_INIT', tC['T_init'][qc] / tU['T_init'][qu], out)
    stats('correction factor C_M', tC['C_M'][qc] / tU['C_M'][qu], out)
    stats('corrected mass M_CORR', tC['M_corr'][qc] / tU['M_corr'][qu], out)
    stats('C_M_3D', tC['C_M_3D'][qc] / tU['C_M_3D'][qu], out)
    out += ['', 'Ratios: 10th, 50th (median), and 90th percentiles over the matched sources.']

    txt = '\n'.join(out)
    print(txt)
    if a.out:
        open(a.out, 'w').write(txt + '\n')


if __name__ == '__main__':
    main()
