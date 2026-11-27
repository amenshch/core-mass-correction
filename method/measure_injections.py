#!/usr/bin/env python3
"""
measure_injections.py -- measure injected cores directly in one or more surface
density maps, at positions that are known, without running an extraction.

For every injected core the script reports, in each map supplied:

  peak      the largest surface density within --peak-radius of the given
            position, in cm^-2;
  bg        the local background, the median surface density in an annulus
            between --bg-inner and --bg-outer, in cm^-2;
  excess    peak minus bg, that is the background-subtracted peak surface
            density in the observational convention;
  mass      the excess integrated over a disk of radius --aperture, converted
            to a gas mass in solar masses at the assumed distance.

The positions are read from a table of injected cores, one line per core, with
the pixel coordinates in two of its columns.  Column numbers are counted from
one, as in the getsf catalogs.  All other columns of that table are copied
through to the output, so that the measurements can be compared directly with
the true masses and contrasts of the injected models.

Example
-------
    python3 measure_injections.py --truth truth_table.txt --x-col 1 --y-col 2 \\
        --map hires=hi.surface.density.r13p5.fits \\
        --map perscale=aquila.surfden.r13p5.fits \\
        --out injections_measured.csv

The pixel coordinates are assumed to follow the FITS convention, the centre of
the first pixel being at 1.0; use --zero-based if they start at 0.
"""

import argparse
import sys

import numpy as np

try:
    from astropy.io import fits
except ImportError:                                              # pragma: no cover
    sys.exit("ERROR: this script requires astropy (pip install astropy)")

PC_CM = 3.0856775814913673e18
AMU = 1.6605402e-24
MU_H2 = 2.8
MSUN = 1.98892e33


def read_truth(path, xcol, ycol, zero_based):
    """Return the pixel coordinates and the raw lines of the table."""
    xs, ys, lines, header = [], [], [], None
    for raw in open(path):
        if raw.lstrip().startswith("#"):
            if header is None:
                header = raw.rstrip()
            continue
        f = raw.split()
        if len(f) < max(xcol, ycol):
            continue
        try:
            x = float(f[xcol - 1])
            y = float(f[ycol - 1])
        except ValueError:
            continue
        xs.append(x if zero_based else x - 1.0)
        ys.append(y if zero_based else y - 1.0)
        lines.append(raw.rstrip())
    return np.array(xs), np.array(ys), lines, header


def pixel_size(header):
    for key in ("CDELT2", "CD2_2", "CDELT1"):
        if key in header:
            return abs(float(header[key])) * 3600.0
    sys.exit("ERROR: no CDELT2 in the header; give --pixel explicitly")


def measure(data, x, y, pix_as, peak_r, bg_in, bg_out, aperture):
    """Measure one source; returns peak, background, excess and aperture sum."""
    n_y, n_x = data.shape
    rad = int(np.ceil(bg_out / pix_as)) + 2
    i0, i1 = int(round(x)) - rad, int(round(x)) + rad + 1
    j0, j1 = int(round(y)) - rad, int(round(y)) + rad + 1
    if i0 < 0 or j0 < 0 or i1 > n_x or j1 > n_y:
        return (np.nan,) * 4
    sub = data[j0:j1, i0:i1]
    yy, xx = np.mgrid[j0:j1, i0:i1]
    r = np.hypot(xx - x, yy - y) * pix_as
    finite = np.isfinite(sub)
    core = finite & (r <= peak_r)
    ann = finite & (r >= bg_in) & (r <= bg_out)
    ap = finite & (r <= aperture)
    if core.sum() == 0 or ann.sum() < 20 or ap.sum() == 0:
        return (np.nan,) * 4
    peak = np.max(sub[core])
    bg = np.median(sub[ann])
    total = float(np.sum(sub[ap] - bg))
    return peak, bg, peak - bg, total


def main():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--truth", required=True,
                   help="table of injected cores")
    p.add_argument("--map", action="append", required=True, metavar="NAME=FILE",
                   help="a surface density map to measure; may be repeated")
    p.add_argument("--out", required=True, help="output comma-separated file")
    p.add_argument("--x-col", type=int, default=1,
                   help="column of the x pixel coordinate, counted from one")
    p.add_argument("--y-col", type=int, default=2,
                   help="column of the y pixel coordinate, counted from one")
    p.add_argument("--zero-based", action="store_true",
                   help="the pixel coordinates in the table start at zero")
    p.add_argument("--pixel", type=float, default=None,
                   help="pixel size in arcsec; taken from the header if absent")
    p.add_argument("--distance", type=float, default=260.0,
                   help="assumed distance in parsec (default: %(default)s)")
    p.add_argument("--peak-radius", type=float, default=9.0,
                   help="radius in arcsec within which the peak is sought")
    p.add_argument("--aperture", type=float, default=20.0,
                   help="radius in arcsec over which the mass is integrated")
    p.add_argument("--bg-inner", type=float, default=60.0,
                   help="inner radius in arcsec of the background annulus")
    p.add_argument("--bg-outer", type=float, default=90.0,
                   help="outer radius in arcsec of the background annulus")
    a = p.parse_args()

    maps = {}
    for entry in a.map:
        if "=" not in entry:
            sys.exit("ERROR: --map expects NAME=FILE, got %s" % entry)
        name, path = entry.split("=", 1)
        with fits.open(path) as hdul:
            hdu = next(h for h in hdul if h.data is not None)
            maps[name] = (np.array(hdu.data, float).squeeze(), hdu.header)
        print("%-12s %s  shape %s" % (name, path, maps[name][0].shape))

    x, y, lines, header = read_truth(a.truth, a.x_col, a.y_col, a.zero_based)
    print("injected cores in the table: %d" % len(x))

    first = next(iter(maps.values()))[1]
    pix = a.pixel if a.pixel else pixel_size(first)
    pix_cm = pix * a.distance * PC_CM / 206264.806
    pix_area = pix_cm ** 2
    to_msun = MU_H2 * AMU * pix_area / MSUN
    print("pixel %.3f arcsec, distance %.1f pc, one pixel of 1e21 cm^-2 "
          "is %.3e solar masses" % (pix, a.distance, 1e21 * to_msun))

    names = sorted(maps)
    cols = []
    for nm in names:
        cols += ["%s_peak" % nm, "%s_bg" % nm, "%s_excess" % nm, "%s_mass" % nm]
    with open(a.out, "w") as fh:
        if header:
            fh.write("# original table header: %s\n" % header)
        fh.write("# x_pix,y_pix," + ",".join(cols) + ",original_line\n")
        nok = 0
        for k in range(len(x)):
            row = [x[k] + 1.0, y[k] + 1.0]
            good = True
            for nm in names:
                data, _ = maps[nm]
                pk, bg, ex, tot = measure(data, x[k], y[k], pix,
                                          a.peak_radius, a.bg_inner,
                                          a.bg_outer, a.aperture)
                row += [pk, bg, ex, tot * to_msun]
                good = good and np.isfinite(pk)
            nok += good
            fh.write(",".join("%g" % v for v in row) + ',"'
                     + lines[k].replace('"', "'") + '"\n')
    print("measured %d of %d cores in every map -> %s" % (nok, len(x), a.out))

    # a short summary, when two maps are supplied
    if len(names) == 2:
        i1 = 2 + cols.index("%s_excess" % names[0])
        i2 = 2 + cols.index("%s_excess" % names[1])
        j1 = 2 + cols.index("%s_mass" % names[0])
        j2 = 2 + cols.index("%s_mass" % names[1])
        vals = []
        for line in open(a.out):
            if line.startswith("#"):
                continue
            f = line.split(",")
            try:
                vals.append([float(f[i1]), float(f[i2]),
                             float(f[j1]), float(f[j2])])
            except ValueError:
                continue
        v = np.array(vals)
        for k, lab in ((0, "background-subtracted peak"), (2, "aperture mass")):
            m = np.isfinite(v[:, k]) & np.isfinite(v[:, k + 1]) \
                & (v[:, k] > 0) & (v[:, k + 1] > 0)
            if m.sum() > 3:
                r = v[m, k + 1] / v[m, k]
                print("%s, %s / %s over %d cores: 16th %.2f, median %.2f, "
                      "84th %.2f" % (lab, names[1], names[0], m.sum(),
                                     *np.percentile(r, [16, 50, 84])))


if __name__ == "__main__":
    main()
