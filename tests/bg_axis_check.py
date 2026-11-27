"""Is the background-column axis of the corrector displaced between training
and application?

On the model side the first observable is SD_emb, the surface density of the
synthetic cloud in which the model was embedded.  On the application side it is
the background surface density measured by getsf on the real hires map beneath
the source.  Two comparisons, both available from files already in hand, test
whether these are the same quantity in practice:

  (a) getsf-measured background against the truth table's local_Sigma, the hires
      value at the placement position.  This isolates the measurement: how well
      does the extraction recover a background it was placed on?

  (b) local_Sigma against the model's nominal SD_emb.  The injections use
      contrast-matched placement, so these need not agree by construction, and
      any systematic difference displaces the query point along the first axis
      relative to where the corrector was trained.

Truth-table columns: 5 local_Sigma, 6 SD_emb, both in cm^-2.
"""
import glob, os, sys
sys.path.insert(0, 'python'); sys.path.insert(0, '.')
import numpy as np
import mass_correction_pipeline_6 as mc
lg = np.log10
DIST, MATCH = 260.0, 10.0
RUNS = ['260811_cmf2p00_n600_s2030', '260812_cmf2p00_n600_s2031',
        '260813_cmf2p00_n600_s2032', '260814_cmf2p00_n600_s2033']
VAR = {'': 'as placed', '_flat': 'flattened', '_iso': 'isolated',
       '_flat_iso': 'flattened and isolated'}
mc.VERBOSE = -1


def q(label, r):
    g = np.isfinite(r) & (r > 0)
    lr = lg(r[g])
    print('     %-46s n=%4d  median %6.3f  tenth %6.3f  ninetieth %6.3f  '
          'typical error %6.3f dex'
          % (label, int(g.sum()), np.median(r[g]),
             *np.percentile(r[g], [10, 90]), np.median(np.abs(lr))),
          flush=True)


for suf, name in VAR.items():
    ms, lsig, sdemb = [], [], []
    for run in RUNS:
        d = os.path.join('Aquila-injections', run + suf)
        cats = [os.path.join(d, 'Aquila.s.sources.ok.cat'),
                os.path.join(d, 'Aquila.s.sources.ok.add.cat')]
        if not os.path.exists(cats[0]):
            continue
        mc.load_getsf._dist = DIST
        g = mc.load_getsf(cats)
        tr = []
        for ln in open(glob.glob(os.path.join(d, '*_truth.txt'))[0]):
            if not ln.strip() or ln.lstrip()[0] in '#!':
                continue
            f = ln.split()
            tr.append((float(f[2]), float(f[3]), float(f[4]), float(f[5])))
        tx = np.array([a[0] for a in tr]); ty = np.array([a[1] for a in tr])
        # load_getsf does not expose pixel coordinates, so they are read from
        # the catalogue directly and matched to the source numbers it returns.
        from getsf_columns import GetsfTable
        t = GetsfTable(cats[0])
        cc = lambda n: np.asarray(t.col(n), float).reshape(t.nrows, -1)[:, 0]
        pos = {int(n): (x, y) for n, x, y in zip(cc('NO'), cc('XCO_P'),
                                                cc('YCO_P'))}
        for i in range(len(g['NO'])):
            if int(g['NO'][i]) not in pos:
                continue
            xi_, yi_ = pos[int(g['NO'][i])]
            dd = np.hypot(tx - xi_, ty - yi_); k = int(np.argmin(dd))
            if dd[k] >= MATCH:
                continue
            ms.append(g['Nbg'][i]); lsig.append(tr[k][2]); sdemb.append(tr[k][3])
        continue
        for i in range(len(g['NO'])):
            dd = np.hypot(tx - xs[i], ty - ys[i]); k = int(np.argmin(dd))
            if dd[k] >= MATCH:
                continue
            ms.append(g['Nbg'][i]); lsig.append(tr[k][2]); sdemb.append(tr[k][3])
    if not ms:
        continue
    ms, lsig, sdemb = map(np.array, (ms, lsig, sdemb))
    print('\n  %s: %d matched sources' % (name, len(ms)), flush=True)
    q('getsf background / hires value at placement', ms / lsig)
    q('hires value at placement / model SD_emb', lsig / sdemb)
    q('getsf background / model SD_emb', ms / sdemb)
