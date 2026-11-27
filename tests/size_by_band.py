"""How the size estimators behave as the source becomes marginally resolved.

The quantity sought throughout is the intrinsic width at half maximum of the
surface-density profile, the axis FWHM_pc of the model converted to arcseconds
at the distance at which the grid was rendered.  What changes between the blocks
below is the measurement from which it is estimated: the width of the truncated,
rim-subtracted profile in one band, at that band's resolution.

    band                 beam O (arcsec)
    surface density      13.5
    250 micrometres      18.2
    350 micrometres      24.9
    500 micrometres      36.3

Three estimators are compared in each band:
    the measured width itself, undeconvolved;
    the Gaussian deconvolution sqrt(H^2 - O^2), standard practice;
    a grid inversion trained on the background surface density, the
      concentration, the footprint factor, the footprint radius and that band's
      measured width, with one node withheld at a time.

Results are also split by the resolvedness H/O, since the interest of a better
estimator lies in the marginally resolved regime where H/O approaches unity and
the Gaussian deconvolution becomes ill-conditioned.
"""
import sys
sys.path.insert(0, 'python')
import numpy as np
import mass_correction_pipeline_6 as mc
lg = np.log10
BEAM = {'SD': 13.5, '250': 18.2, '350': 24.9, '500': 36.3}
KEYW = {'SD': 'w_sd', '250': 'w250', '350': 'w350', '500': 'w500'}
HO_BINS = [1.0, 1.3, 1.6, 2.0, 3.0, 10.0]


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


def line(label, ratio, ho):
    g = np.isfinite(ratio) & (ratio > 0)
    out = ('     %-32s n=%5d  median %6.3f  typical error %6.3f dex  |'
           % (label, int(g.sum()), np.median(ratio[g]),
              np.median(np.abs(lg(ratio[g])))))
    for lo, hi in zip(HO_BINS[:-1], HO_BINS[1:]):
        m = g & (ho >= lo) & (ho < hi)
        out += (' %6.3f' % np.median(np.abs(lg(ratio[m])))
                if m.sum() > 30 else '      -')
    print(out, flush=True)


mc.VERBOSE = 0
S = mc.build_recovery_samples('cats/bes_model_grid_final3_catalog',
                              'cats/bes_model_grid_final3_recovery_tables')
ok = np.isfinite(S['rbe_pc']) & (S['rbe_pc'] > 0)
node = np.asarray(S['node'])[ok]
key = node[:, 0] * 100000 + node[:, 1] * 1000 + node[:, 2]
G3 = read_grid('cats/bes_grid_final3.txt')
SC = 206264.806 / mc.GRID_DISTANCE_PC
truth = np.array([G3[tuple(n)]['FWHM_pc'] * SC if tuple(n) in G3 else np.nan
                  for n in node])
sig, conc, phi = S['sig'][ok], S['conc'][ok], S['cfoot'][ok]
r_as = phi * S['w_sd'][ok] / 2.0

print('Typical error, in decimal exponents, of each estimate of the intrinsic '
      'half-maximum\nwidth.  Columns after the bar give it in bins of the '
      'resolvedness H/O = %s.\n'
      % ', '.join('%.1f-%.1f' % (a, b)
                  for a, b in zip(HO_BINS[:-1], HO_BINS[1:])), flush=True)

for b in ('SD', '250', '350', '500'):
    w = S[KEYW[b]][ok]
    X = np.column_stack([lg(sig), lg(conc), lg(phi), lg(r_as), lg(w)])
    g = (np.all(np.isfinite(X), axis=1) & np.isfinite(truth) & (truth > 0)
         & (w > 0))
    Xb, kb, tb, wb = X[g], key[g], truth[g], w[g]
    ho = wb / BEAM[b]
    pred = np.full(len(kb), np.nan)
    for k in np.unique(kb):
        o = kb == k
        keep = ~o
        if keep.sum() < 100:
            continue
        pred[o] = mc.LocalLinear(Xb[keep], lg(tb[keep]))(Xb[o])
    pred = 10 ** pred
    dec = np.sqrt(np.clip(wb ** 2 - BEAM[b] ** 2, 0, None))
    print('  %s band, beam %.1f arcsec: %d samples, median H/O %.2f, '
          'unresolved (H/O<1) %.0f%%'
          % (b, BEAM[b], len(wb), np.median(ho), 100 * np.mean(ho < 1)),
          flush=True)
    line('measured width, undeconvolved', wb / tb, ho)
    line('Gaussian deconvolution', dec / tb, ho)
    line('grid inversion, node withheld', pred / tb, ho)
    print(flush=True)
