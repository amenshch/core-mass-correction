"""Separating what the low-resolution width contributes from what the
surface-density observables contribute.

In the band-by-band test the inversion received Sigma, zeta, phi and the
footprint radius from the surface-density band, and only the width from the
band under test.  It therefore measured the realistic Herschel case, in which a
surface-density map exists alongside the photometry, but it could not show how
much of the accuracy came from the low-resolution width itself.

Here the inversion is given only two inputs: the background surface density,
which any map supplies, and the width measured in the band under test.  The
target is unchanged, the intrinsic width at half maximum of the surface-density
profile.  One node is withheld at a time.  Comparing these figures with the
five-input ones separates the two contributions.
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
    out = ('     %-34s n=%5d  median %6.3f  typical error %6.3f dex  |'
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
sig = S['sig'][ok]

print('Two-input inversion: background surface density and the width measured '
      'in the band\nunder test, nothing else.  Columns after the bar give the '
      'typical error in bins of\nthe resolvedness H/O = %s.\n'
      % ', '.join('%.1f-%.1f' % (a, b)
                  for a, b in zip(HO_BINS[:-1], HO_BINS[1:])), flush=True)

for b in ('SD', '250', '350', '500'):
    w = S[KEYW[b]][ok]
    X = np.column_stack([lg(sig), lg(w)])
    g = (np.all(np.isfinite(X), axis=1) & np.isfinite(truth) & (truth > 0)
         & (w > 0))
    Xb, kb, tb, wb = X[g], key[g], truth[g], w[g]
    ho = wb / BEAM[b]
    pred = np.full(len(kb), np.nan)
    for k in np.unique(kb):
        o = kb == k
        if (~o).sum() < 100:
            continue
        pred[o] = mc.LocalLinear(Xb[~o], lg(tb[~o]))(Xb[o])
    pred = 10 ** pred
    print('  %s band, beam %.1f arcsec: %d samples' % (b, BEAM[b], len(wb)),
          flush=True)
    line('two-input inversion', pred / tb, ho)
