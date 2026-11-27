"""Why the temperature factor C_T has no effect on injected sources.

Compares, per source, the correction the model grid predicts with the
correction the source actually needs, the latter being the ratio of the
injected model's mass-averaged dust temperature Td_avg to the temperature
returned by the fit to the spectral energy distribution of the extracted
source.  If the predicted factors cluster near unity while the required ones
scatter widely, the corrector is declaring that no correction is needed and the
real-data error has a different origin from the one the models contain.
"""
import glob, os, sys
sys.path.insert(0, 'python'); sys.path.insert(0, '.')
import numpy as np
import mass_correction_pipeline_6 as mc
from getsf_columns import GetsfTable
from scipy.stats import spearmanr
lg = np.log10
DIST, MATCH = 260.0, 10.0
RUNS = ['260811_cmf2p00_n600_s2030', '260812_cmf2p00_n600_s2031',
        '260813_cmf2p00_n600_s2032', '260814_cmf2p00_n600_s2033']
VAR = {'': 'as placed', '_flat': 'flattened'}


def cols(path, c):
    keys, out = [], {n: [] for n in c}
    for l in open(path):
        if l.startswith('#') or not l.strip():
            continue
        f = l.split(); keys.append((int(f[1]), int(f[2]), int(f[3])))
        for n, i in c.items():
            try: out[n].append(float(f[i - 1]))
            except Exception: out[n].append(np.nan)
    return {k: {n: out[n][i] for n in c} for i, k in enumerate(keys)}


mc.VERBOSE = 0
S = mc.build_recovery_samples('cats/bes_model_grid_final3_catalog',
                              'cats/bes_model_grid_final3_recovery_tables')
ok = np.isfinite(S['rbe_pc']) & (S['rbe_pc'] > 0)
th = S['rbe_pc'][ok] * 206264.806 / mc.GRID_DISTANCE_PC
X = np.column_stack([lg(S['sig'][ok]), lg(S['conc'][ok]), lg(S['cfoot'][ok]),
                     lg(th)])
node = np.asarray(S['node'])[ok]
c3 = cols('cats/bes_model_grid_final3_catalog', {'Td': 12, 'T4': 17, 'T3': 23})
g = lambda n: np.array([c3[tuple(x)][n] if tuple(x) in c3 else np.nan
                        for x in node])
ts = np.where(np.isfinite(g('T4')) & (g('T4') > 0), g('T4'), g('T3'))
ct = g('Td') / ts
m = np.all(np.isfinite(X), axis=1) & np.isfinite(ct) & (ct > 0)
print('model-side C_T = Td_avg/T_SED over %d samples: median %.4f, '
      'tenth to ninetieth percentile %.4f to %.4f'
      % (m.sum(), np.median(ct[m]), *np.percentile(ct[m], [10, 90])),
      flush=True)
est = mc.LocalLinear(X[m], lg(ct[m]))

c2 = cols('cats/bes_model_grid_final2_catalog', {'Td': 12})
for suf, name in VAR.items():
    rows = []
    for r in RUNS:
        d = os.path.join('Aquila-injections', r + suf)
        cats = [os.path.join(d, 'Aquila.s.sources.ok.cat'),
                os.path.join(d, 'Aquila.s.sources.ok.add.cat')]
        if not os.path.exists(cats[0]):
            continue
        tab, _ = mc.correct_cloud(cats, DIST)
        t = GetsfTable(cats[0])
        cc = lambda n: np.asarray(t.col(n), float).reshape(t.nrows, -1)[:, 0]
        idx = {int(n): i for i, n in enumerate(cc('NO'))}
        xs, ys = cc('XCO_P'), cc('YCO_P')
        tr = []
        for ln in open(glob.glob(os.path.join(d, '*_truth.txt'))[0]):
            if not ln.strip() or ln.lstrip()[0] in '#!':
                continue
            f = ln.split(); tg = f[1]
            tr.append((float(f[2]), float(f[3]),
                       (int(tg[1:3]), int(tg[4:6]), int(tg[7:9]))))
        tx = np.array([q[0] for q in tr]); ty = np.array([q[1] for q in tr])
        for i, no in enumerate(tab['NO']):
            j = idx.get(int(no))
            if j is None:
                continue
            dd = np.hypot(tx - xs[j], ty - ys[j]); k = int(np.argmin(dd))
            if dd[k] >= MATCH or tr[k][2] not in c2:
                continue
            rows.append((tab['sigma'][i], tab['conc'][i], tab['rfoot'][i],
                         tab['foot_as'][i], tab['T_init'][i],
                         c2[tr[k][2]]['Td']))
    D = np.array(rows)
    q = np.column_stack([lg(D[:, 0]), lg(D[:, 1]), lg(D[:, 2]), lg(D[:, 3])])
    v = np.all(np.isfinite(q), axis=1) & np.isfinite(D[:, 4]) & (D[:, 4] > 0)
    pred = 10 ** est(q[v])
    need = D[v, 5] / D[v, 4]
    print('\n%s: %d sources' % (name, v.sum()))
    print('  predicted C_T : median %.4f, tenth to ninetieth %.4f to %.4f'
          % (np.median(pred), *np.percentile(pred, [10, 90])))
    print('  required C_T  : median %.4f, tenth to ninetieth %.4f to %.4f'
          % (np.median(need), *np.percentile(need, [10, 90])))
    print('  rank correlation of predicted with required: %+.3f'
          % spearmanr(pred, need).statistic)
