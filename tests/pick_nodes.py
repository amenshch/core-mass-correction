"""Choose a small set of final3 grid nodes at which to compare the hires-derived
surface-density profile with the model's true one.

The exposure of the corrector to hires lies in the radial gradient of the bias
inside a core, so the set must span the quantities on which that gradient can
depend: the embedding background surface density, which sets how strongly the
core stands out; the truncation parameter xi_max, which sets the concentration
and hence the steepness of the internal temperature gradient; and the observed
angular size, which sets how many resolution elements the profile is sampled
with.  Nodes typical of the real Aquila population are included so that the
diagnostic covers the regime most catalogue sources occupy, and hole nodes,
whose outer density falls below the cloud, because those are the ones the fourth
observable exists to rescue.
"""
import numpy as np
GRID = 'cats/bes_grid_final3.txt'
CAT = 'cats/bes_model_grid_final3_catalog'
XI_CRIT = 6.451
COLS = {'M_BE': 8, 'FWHMbs': 54, 'FWHMbsl': 61}


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


def read_cat(p, cols):
    out = {}
    for l in open(p):
        if l.startswith('#') or not l.strip():
            continue
        f = l.split()
        d = {}
        for n, i in cols.items():
            try: d[n] = float(f[i - 1])
            except Exception: d[n] = np.nan
        out[(int(f[1]), int(f[2]), int(f[3]))] = d
    return out


G, C = read_grid(GRID), read_cat(CAT, COLS)
keys = [k for k in G if k in C and G[k]['detect'] == 1
        and np.isfinite(C[k]['FWHMbsl']) and C[k]['FWHMbsl'] > 0]
sd = np.array([G[k]['SD_emb'] for k in keys])
xi = np.array([G[k]['xi_max'] for k in keys])
con = np.array([G[k]['contrast'] for k in keys])
hol = np.array([G[k]['hole'] for k in keys])
size = np.array([C[k]['FWHMbsl'] for k in keys])
mbe = np.array([G[k]['M_BE'] for k in keys])
sdvals = np.unique(sd)
print('detectable nodes with a finite observed size: %d' % len(keys))
print('background column values available: %s\n'
      % ' '.join('%.3g' % v for v in sdvals))

chosen, why = [], {}


def take(mask, label, n=1):
    idx = np.where(mask)[0]
    if len(idx) == 0:
        return
    for i in idx[:n]:
        k = keys[i]
        if k not in chosen:
            chosen.append(k); why[k] = label


# 1. three background columns x three stability classes, at mid observed size
smid = (size > 25) & (size < 60)
for sdl, sdlab in ((sdvals[1], 'low'), (sdvals[4], 'mid'), (sdvals[7], 'high')):
    m0 = smid & (np.abs(sd / sdl - 1) < 1e-6)
    for lab, mm in (('sub-critical', xi < 0.7 * XI_CRIT),
                    ('near-critical', np.abs(xi / XI_CRIT - 1) < 0.15),
                    ('super-critical', xi > 3 * XI_CRIT)):
        take(m0 & mm, '%s background, %s' % (sdlab, lab))

# 2. size ladder at mid background, near-critical
m0 = (np.abs(sd / sdvals[4] - 1) < 1e-6) & (np.abs(xi / XI_CRIT - 1) < 0.25)
for lo, hi, lab in ((10, 20, 'compact'), (35, 50, 'medium'),
                    (80, 200, 'extended')):
    take(m0 & (size > lo) & (size < hi), 'mid background, %s size' % lab)

# 3. nodes typical of the real Aquila population
d = np.hypot((np.log10(con / 1.235)) / 0.05, (np.log10(size / 40.0)) / 0.10)
for i in np.argsort(d)[:2]:
    k = keys[i]
    if k not in chosen:
        chosen.append(k); why[k] = 'typical of real Aquila sources'

# 4. hole nodes
take((hol == 1) & smid & (xi > XI_CRIT), 'hole node, super-critical')
take((hol == 1) & smid & (xi < XI_CRIT), 'hole node, sub-critical')

print('%-14s %11s %8s %9s %8s %8s %6s  %s'
      % ('node', 'SD_emb', 'M_BE', 'xi_max', 'size', 'contrast', 'hole',
         'reason for inclusion'))
for k in chosen:
    print('i%02dj%02dk%02d %11.3g %8.4g %9.3f %8.1f %8.3f %6d  %s'
          % (k[0], k[1], k[2], G[k]['SD_emb'], G[k]['M_BE'], G[k]['xi_max'],
             C[k]['FWHMbsl'], G[k]['contrast'], G[k]['hole'], why[k]))
print('\nunits: SD_emb in cm^-2, M_BE in solar masses, size = FWHMSDbsl in '
      'arcsec,\ncontrast = peak/background on the surface-density band; '
      'xi_max critical at 6.451')
