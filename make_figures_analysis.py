#!/usr/bin/env python3
# make_figures_analysis.py - the remaining document figures (Figs 2, 4, 6, 8, 9, 10, 11, 12, 14, 16, 17).
# Companion of make_figures.py (re-uses its Fig / panel / legend helpers, so the style and the layout checks are identical).
#
# EVERY plotted number is computed by the package's own programs, nothing is typed in:
#   ./rpt_adapt            -> all deterministic runs (through analysis_scripts.run / cT_star, as in Sections 10-12)
#   analysis_scripts.mc    -> Monte Carlo establishment probability (Fig 16a)       analysis_scripts.q_theory -> Eq 26 (Fig 16)
#   analysis_scripts.tail_sim -> Gillespie tail simulation (Fig 17)
#   screen_report.txt      -> the screened-best profile written by ./engine --screen (Fig 4)
# Each figure's profiles, grids and parameter lists are defined in the function that draws it and are stated in the figure's caption.
# figure_data/<name>.csv holds the plotted numbers.
#
# Usage:  g++ -std=c++20 -O2 rpt_adapt.cpp -o rpt_adapt ; ./engine --screen screen.txt input.txt errors.txt   (once)
#         python3 make_figures_analysis.py [outdir] [fig2 fig4 fig6 fig8 fig9 fig10 fig11 fig12 fig14 fig16 fig17]   (default: all)
import math, os, sys
import numpy as np
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import make_figures as mf
AS, Fig, panel, legend, rich = mf.AS, mf.Fig, mf.panel, mf.legend, mf.rich
OUT = mf.OUT
INP = {k: v for k, v in mf.read_input(os.path.join(HERE, 'input.txt')).items() if k not in ('w1', 'w2', 'w3', 'w4', 'C_crit', 'H_crit')}
BEST = dict(AS.BEST)                                         # screened-best sampled profile (analysis_scripts.py)
BLUE, CORAL, GREEN, GREY, PURPLE, GOLD, RED = '#185FA5', '#D85A30', '#1D9E75', '#6B6B6B', '#7F4FA3', '#B8860B', '#A32D2D'
CW_BAND = (5.0, 15.0)                                        # WT clearance range c_W (document, Section 5)
BAND_COL = '#8a8a8a'

def run(o, adaptive=0): return AS.run({**o, 'adaptive': adaptive})
def tstar(prof, extra=None, adaptive=0, lo=0.003, hi=20.0):
    r = AS.cT_star(prof, extra, lo=lo, hi=hi, iters=22, adaptive=adaptive)
    return float('nan') if (r != r or r >= hi * 0.999) else r      # no extinction down to lo, or censored at hi -> nan
def pmap(f, items):
    with Pool(2) as p: return p.map(f, items, chunksize=4)
def segments(x, y, color, dash=None, width=2.4):
    """Split a curve with gaps (nan) into polylines; isolated points get a dot marker."""
    out, cur = [], []
    for a, b in list(zip(x, y)) + [(None, float('nan'))]:
        if b == b: cur.append((a, b))
        else:
            if len(cur) > 1: out.append(dict(x=[c[0] for c in cur], y=[c[1] for c in cur], color=color, dash=dash, width=width))
            elif len(cur) == 1: out.append(dict(x=[cur[0][0]], y=[cur[0][1]], color=color, noline=True, marker='dot'))
            cur = []
    return out
def vband(F, X, Y, W, H, px, lo, hi, R=70):
    L, T, B = 96, 36, 56
    F.add(f"<rect x='{px(lo):.1f}' y='{Y + T}' width='{px(hi) - px(lo):.1f}' height='{H - T - B}' fill='{BAND_COL}' fill-opacity='0.28'/>")
def finish(F, name):
    n = F.check(); F.save(os.path.join(OUT, name)); return n
CT_LAB = 'free-TIP clearance *c*_{T} (d^{−1})'

# ------------------------------------------------------------------------------------------------ simple bar helpers
def axis_range(vals, log, extra=()):
    vals = [v for v in list(vals) + list(extra) if v == v]
    if log:
        lo, hi = math.floor(math.log10(min(vals)) - 0.3), math.ceil(math.log10(max(vals)) - 1e-9)
        return lo, (hi if hi > lo else lo + 1), 1
    st = mf.nice_step(0, max(vals) * 1.06, 5); return 0.0, math.ceil(max(vals) * 1.06 / st - 1e-9) * st, st
def bars(F, X, Y, W, H, title, ylab, cats, vals, cols, log=False, hline=None, hlab=None):
    L, R, T, B = 96, 70, 36, 56; pw, ph = W - L - R, H - T - B
    lo, hi, st = axis_range(vals, log, [hline] if hline else [])
    yy = lambda v: Y + T + (1 - ((math.log10(max(v, 1e-300)) if log else v) - lo) / (hi - lo)) * ph
    F.text(X + L, Y + 22, title, mf.S_PTITLE, bold=True, fill='#1a1a1a', kind='title')
    F.add(f"<rect x='{X + L}' y='{Y + T}' width='{pw}' height='{ph}' fill='#fafafa'/>")
    v = math.ceil(lo / st - 1e-9) * st
    while v <= hi + st * 1e-9:
        yp = Y + T + (1 - (v - lo) / (hi - lo)) * ph
        F.line(X + L, yp, X + L + pw, yp, '#e4e4e4', 1, check=False); F.line(X + L - 5, yp, X + L, yp, '#444', 1.2, check=False)
        if log:
            m, pl, w = mf.pow10(int(round(v))); F.add(f"<text x='{X + L - 9:.1f}' y='{yp + 4.5:.1f}' font-size='13' text-anchor='end' fill='#333'>{m}</text>")
            F.texts.append(dict(bb=(X + L - 9 - w, yp + 4.5 - 10.4, X + L - 9, yp + 8.5), plain=pl, kind='tick', note=False, owner=None, w=w))
        else: F.text(X + L - 9, yp + 4.5, mf.fmt_tick(v, st), 13, 'end', fill='#333', kind='tick')
        v += st
    slot = pw / len(vals)
    for i, (c, val, col) in enumerate(zip(cats, vals, cols)):
        x0 = X + L + slot * (i + 0.2); top = yy(val); base = Y + T + ph
        F.add(f"<rect x='{x0:.1f}' y='{top:.1f}' width='{slot * 0.6:.1f}' height='{max(base - top, 0.5):.1f}' fill='{col}'/>")
        F.text(X + L + slot * (i + 0.5), base + 21, c, 13, 'middle', fill='#333', kind='tick')
    F.line(X + L, Y + T + ph, X + L + pw, Y + T + ph, '#222', 1.6, check=False); F.line(X + L, Y + T, X + L, Y + T + ph, '#222', 1.6, check=False)
    F.text(X + 18, Y + T + ph / 2 + 5, ylab, mf.S_LABEL, 'middle', fill='#333', kind='axis', rot=-90)
    if hline:
        F.line(X + L, yy(hline), X + L + pw, yy(hline), RED, 1.5, '2 4', check=False); F.text(X + L + pw + 8, yy(hline) + 4.5, hlab, 13, fill=RED, kind='ref', note=True)

# ================================================================================================ Fig 6
def fig6():
    base = {**INP, 'V_W0': 1e-12, 'tip_dose': 0, 'ext_floor': 0, 't_end': 60, 'dt_out': 0.25}
    Ts = INP['lambda'] / INP['d_T']; Cs = INP['C_rest']; cyt = INP['delta_N'] * Cs / (INP['K_F'] + Cs)
    def eig(eAb, eCl):                                       # linearised (E_W, I_W, V_W) block at the infection-free state, from the right-hand sides of rpt_adapt.cpp
        bt = (1 - eAb) * INP['beta_W'] * Ts; k = INP['k_W']; d = INP['delta_W'] + cyt; p = INP['p_W']; c = INP['c_W'] * (1 + INP['kappa'] * eCl)
        return max(np.linalg.eigvals(np.array([[-k, 0, bt], [k, -d, 0], [0, p, -c]])).real)
    R0 = run({**base})['R0']
    def sim(over):                                           # growth rate of V_W in the linear regime (t >= 3 d, V_W < 1e4; V_W(0) = 1e-12 so no target cells are depleted)
        _, c = mf.rpt_traj({**base, **over}); t, V = np.array(c['t']), np.array(c['VW']); m = (t >= 3) & (V < 1e4) & (V > 0)
        return float(np.polyfit(t[m] - t[m].mean(), np.log(V[m]), 1)[0])
    ga, gb = np.linspace(0, 0.98, 25), np.linspace(0, 12, 25)
    ra = pmap(sim_a, [(base, e) for e in ga]); rb = pmap(sim_b, [(base, e) for e in gb])
    fa, fb = np.linspace(0, 0.99, 200), np.linspace(0, 12, 200)
    ea, eb = [eig(e, 0) for e in fa], [eig(0, e) for e in fb]
    err = max(max(abs(r - eig(e, 0)) for e, r in zip(ga, ra)), max(abs(r - eig(0, e)) for e, r in zip(gb, rb)))
    mf.write_csv('fig06_growth_rate', ['panel', 'effect', 'sim_rate', 'eigenvalue'], [('a_entry_blocking', e, r, eig(e, 0)) for e, r in zip(ga, ra)] + [('b_clearance_gain', e, r, eig(0, e)) for e, r in zip(gb, rb)])
    F = Fig('fig06', 1260, 540); F.text(14, 30, 'Threshold test at the effective reproduction number R_{eff} = 1', mf.S_TITLE, bold=True, kind='title')
    for k, (ttl, xl, xs, ys, fx, fy, thr, lab, xhi) in enumerate([('(a) Entry blocking', 'entry-blocking effect *ε*_{Ab}', ga, ra, fa, ea, 1 - 1 / R0, '1 − 1/*R*_{0}', 1.0),
                                                                    ('(b) Clearance gain', 'clearance gain *ε*_{Clr}', gb, rb, fb, eb, R0 - 1, '*R*_{0} − 1', 12.0)]):
        px, py = panel(F, 8 + 628 * k, 46, 620, 420, ttl, dict(label=xl, lo=0, hi=xhi), dict(label='growth rate of *V*_{W} (d^{−1})'),
                       [dict(x=list(fx), y=fy, color=BLUE), dict(x=list(xs), y=ys, color=CORAL, noline=True, marker='dot')],
                       hlines=[(0.0, '0', '#444', '2 4')], vlines=[(thr, lab, RED, '6 4')])
    legend(F, 14, 500, [('largest eigenvalue of Eq 16', BLUE, None), ('full nonlinear simulation', CORAL, None, 'dot'), ('analytic threshold', RED, '6 4')], 1200)
    print(f'fig06: max |simulated - eigenvalue| growth rate = {err:.2e} d^-1 (document: 8e-5); thresholds 1-1/R0 = {1 - 1 / R0:.5f}, R0-1 = {R0 - 1:.4f} (document 0.88545, 7.7302)')
    return finish(F, 'fig06_threshold_growth_rate.svg')
def sim_a(a): base, e = a; return _sim(base, {'eps_Ab': e})
def sim_b(a): base, e = a; return _sim(base, {'eps_Clr': e})
def _sim(base, over):
    _, c = mf.rpt_traj({**base, **over}); t, V = np.array(c['t']), np.array(c['VW']); m = (t >= 3) & (V < 1e4) & (V > 0)
    return float(np.polyfit(t[m] - t[m].mean(), np.log(V[m]), 1)[0])

# ================================================================================================ Fig 2
def _ext_time(a):
    prof, ct, primed = a; s = run({**prof, 'c_T': ct, 'primed': primed}); return s['tExt'] if s['extinct'] == 1 else float('nan')
WEAK = dict(beta_T=5e-7, alpha=0.2, p_T=1500, tip_dose=1e8, tip_time=1.0)
def profiles():                                              # Section 7.2: strongest / screened best / weaker (other parameters of the screened best)
    # strongest: Nelder-Mead refined profile (analysis_scripts.REFINED) at dose 1e9, day 1 (thresholds 0.46 / 1.15 d^-1)
    # weaker: beta_T 5e-7, alpha 0.2, p_T 1500 at dose 1e8, day 1 (no primed-cell extinction; abortive threshold 0.525 d^-1)
    return [('strongest (dose 10^{9}, day 1)', {**AS.REFINED, 'tip_dose': 1e9, 'tip_time': 1.0}, BLUE), ('screened best', BEST, GREEN), ('weaker (dose 10^{8}, day 1)', WEAK, CORAL)]
def fig2():
    grid = list(np.logspace(math.log10(0.01), math.log10(15), 45)); ser, rows, leg, tri = [], [], [], []
    for name, prof, col in profiles():
        for primed, dash in ((1, None), (0, '8 5')):
            y = pmap(_ext_time, [(prof, c, primed) for c in grid]); rows += [(name, primed, c, v) for c, v in zip(grid, y)]
            ser += segments(grid, y, col, dash)
            ts = tstar(prof, {'primed': primed}); print(f"fig02 {name:32s} {'primed' if primed else 'abortive':8s} c_T* = {ts:.3f}" if ts == ts else f"fig02 {name:32s} {'primed' if primed else 'abortive':8s} c_T*: no extinction down to 0.003")
            if ts == ts: tri.append(dict(x=[ts], y=[1.4], color=col, noline=True, marker='tri'))
        leg.append((name, col, None))
    mf.write_csv('fig02_time_to_extinction', ['profile', 'primed', 'c_T', 'time_to_extinction_days'], rows)
    F = Fig('fig02', 880, 600); F.text(14, 30, 'Time to WT extinction against free-TIP clearance', mf.S_TITLE, bold=True, kind='title')
    px, py = panel(F, 8, 46, 860, 450, '(a) Three TIP profiles', dict(label=CT_LAB, log=True, lo=-2, hi=2), dict(label='time to WT extinction (days)', lo=0, hi=30), ser + tri)
    vband(F, 8, 46, 860, 450, px, *CW_BAND)
    legend(F, 14, 530, leg + [('primed cells tracked (solid)', '#333', None), ('abortive-entry variant (dashed)', '#333', '8 5'), ('reference threshold', '#333', None, 'tri'), ('WT clearance range 5–15 d^{−1}', BAND_COL, None, 'patch')], 840)
    return finish(F, 'fig02_extinction_time_vs_cT.svg')

# ================================================================================================ Fig 4
def screened_best():
    rep = open(os.path.join(HERE, 'screen_report.txt'), encoding='utf8').read().splitlines()
    i = [k for k, l in enumerate(rep) if l.startswith('== Top')][0]; hdr, row = rep[i + 1].split(), rep[i + 2].split()
    names = hdr[2:hdr.index('peakV')]; return dict(zip(names, map(float, row[3:3 + len(names)])))
def fig4():
    sb = screened_best(); print('fig04 screened best from screen_report.txt:', sb, '| analysis_scripts BEST:', BEST)
    arms = [('untreated', {**sb, 'tip_dose': 0}, GREY, '8 5'), ('weak profile', {**sb, **WEAK}, CORAL, None), ('best screened profile', sb, BLUE, None)]
    res = [(n, col, dash) + mf.rpt_traj(o) for n, o, col, dash in arms]
    mf.write_csv('fig04_arms', ['arm', 'peak_V_W', 'peak_C', 'extinct'], [(n, s['peakV'], s['peakC'], s['extinct']) for n, _, _, s, _ in res])
    F = Fig('fig04', 1290, 470); F.text(14, 30, 'Untreated infection against a weak and the best screened TIP profile', mf.S_TITLE, bold=True, kind='title')
    panel(F, 8, 46, 470, 380, '(a) Free WT virus over time', dict(label='time (days)'), dict(label='copies/mL', log=True), [dict(x=c['t'], y=c['VW'], color=col, dash=dash) for _, col, dash, _, c in res],
          hlines=[(2e-4, '*V*_{ext}', RED, '2 4')], R=110)
    cats = ['untreated', 'weak', 'best']; cols = [GREY, CORAL, BLUE]
    bars(F, 484, 46, 400, 380, '(b) Peak WT load', 'peak *V*_{W} (copies/mL)', cats, [s['peakV'] for _, _, _, s, _ in res], cols, log=True)
    bars(F, 884, 46, 400, 380, '(c) Peak cytokine', 'peak cytokine (pg/mL)', cats, [s['peakC'] for _, _, _, s, _ in res], cols, hline=500, hlab='*C*_{crit}')
    legend(F, 14, 450, [('untreated', GREY, '8 5'), ('weak profile (dose 10^{8}, day 1)', CORAL, None), ('best screened profile', BLUE, None)], 1200)
    for n, _, _, s, _ in res: print(f"fig04 {n:22s} peakV={s['peakV']:.3g} peakC={s['peakC']:.3g} extinct={int(s['extinct'])}")
    return finish(F, 'fig04_weak_vs_best.svg')

# ================================================================================================ Fig 8
def _ext(a):
    prof, ct, dose, primed = a; return run({**prof, 'c_T': ct, 'tip_dose': dose, 'primed': primed})['extinct']
def fig8():
    cg = list(np.logspace(math.log10(0.01), math.log10(20), 36)); dg = list(np.logspace(5, 9, 48))       # 36 x 48 = 1728 combinations
    Z = {}
    for primed in (1, 0):
        r = pmap(_ext, [(BEST, c, d, primed) for d in dg for c in cg]); Z[primed] = np.array(r).reshape(len(dg), len(cg))
    mf.write_csv('fig08_extinction_map', ['primed', 'c_T', 'dose', 'extinct'], [(p, c, d, int(Z[p][i, j])) for p in (1, 0) for i, d in enumerate(dg) for j, c in enumerate(cg)])
    F = Fig('fig08', 1290, 640); F.text(14, 30, 'Extinction map in (c_{T}, dose) space for the fixed profile', mf.S_TITLE, bold=True, kind='title')
    lc, ld = np.log10(cg), np.log10(dg)
    for k, (primed, ttl) in enumerate([(1, '(a) Primed cells tracked (reference)'), (0, '(b) Abortive-entry variant')]):
        X, Y, W, H = 8 + 640 * k, 46, 630, 480
        px, py = panel(F, X, Y, W, H, ttl, dict(label=CT_LAB, log=True, lo=-2, hi=2), dict(label='TIP dose (copies/mL)', log=True), [dict(x=[0.01, 20], y=[1e5, 1e9], color='none', width=0)])
        ex = lambda a: [lc[0] - (lc[1] - lc[0]) / 2] + [(a[i] + a[i + 1]) / 2 for i in range(len(a) - 1)] + [a[-1] + (a[-1] - a[-2]) / 2]
        xe, ye = ex(lc), ex(ld)
        for i in range(len(dg)):
            for j in range(len(cg)):
                x0, x1 = px(10 ** max(xe[j], -2)), px(10 ** min(xe[j + 1], 2)); y0, y1 = py(10 ** min(ye[i + 1], 9)), py(10 ** max(ye[i], 5))
                F.add(f"<rect x='{x0:.1f}' y='{y0:.1f}' width='{x1 - x0 + 0.4:.1f}' height='{y1 - y0 + 0.4:.1f}' fill='{'#6BBF8E' if Z[primed][i, j] == 1 else '#EFE3C2'}'/>")
        vband(F, X, Y, W, H, px, *CW_BAND)
        ok = [(cg[j], dg[i]) for i in range(len(dg)) for j in range(len(cg)) if Z[primed][i, j] == 1]
        print(f"fig08 {'primed' if primed else 'abortive'}: {len(ok)} of {Z[primed].size} combinations extinct; largest extinct c_T = {max(c for c, _ in ok):.3g} d^-1")
        for i in (0, len(dg) - 1):
            row = [cg[j] for j in range(len(cg)) if Z[primed][i, j] == 1]; print(f"   dose {dg[i]:.2g}: largest extinct c_T = {max(row) if row else float('nan'):.3g}")
    legend(F, 14, 560, [('WT extinct within 30 days', '#6BBF8E', None, 'patch'), ('no extinction', '#EFE3C2', None, 'patch'), ('WT clearance range 5–15 d^{−1}', BAND_COL, None, 'patch')], 1200)
    return finish(F, 'fig08_extinction_map.svg')

# ================================================================================================ Fig 9
def fig9():
    cts = [0.05, 0.1, 0.3, 0.5, 1.0, 2.0, 5.0]                                                         # seven values straddling c_T* = 0.256: 0.05 and 0.1 go extinct, larger values regrow
    cols = ['#08519c', '#3182bd', '#41ab5d', '#d9a400', '#f16913', '#cb181d', '#67000d']
    res = [mf.rpt_traj({**BEST, 'c_T': c, 'dt_out': 0.1}) for c in cts]
    mf.write_csv('fig09_trajectories', ['c_T', 'time_days', 'V_W', 'V_T'], [(c, t, v, w) for c, (_, tr) in zip(cts, res) for t, v, w in zip(tr['t'], tr['VW'], tr['VT'])])
    F = Fig('fig09', 1290, 860); F.text(14, 30, 'Trajectories for seven free-TIP clearance rates', mf.S_TITLE, bold=True, kind='title')
    panel(F, 8, 46, 640, 380, '(a) Free WT virus', dict(label='time (days)'), dict(label='copies/mL', log=True), [dict(x=tr['t'], y=tr['VW'], color=c) for c, (_, tr) in zip(cols, res)], hlines=[(2e-4, '*V*_{ext}', RED, '2 4')])
    panel(F, 650, 46, 640, 380, '(b) Free TIP', dict(label='time (days)'), dict(label='copies/mL', log=True), [dict(x=tr['t'], y=tr['VT'], color=c) for c, (_, tr) in zip(cols, res)])
    panel(F, 8, 436, 640, 380, '(c) Phase portrait', dict(label='free TIP *V*_{T} (copies/mL)', log=True), dict(label='free WT virus *V*_{W} (copies/mL)', log=True), [dict(x=tr['VT'], y=tr['VW'], color=c) for c, (_, tr) in zip(cols, res)])
    legend(F, 700, 500, [(f'*c*_{{T}} = {c:g} d^{{−1}}', col, None) for c, col in zip(cts, cols)], 560)
    for c, (s, tr) in zip(cts, res):
        v = np.array(tr['VW']); t = np.array(tr['t']); pos = v > 0
        print(f"fig09 c_T={c:<5g} extinct={int(s['extinct'])} tExt={s['tExt']:.3g} min V_W after dose={v[(t > 0.5) & pos].min():.3g}")
    return finish(F, 'fig09_trajectories_near_threshold.svg')

# ================================================================================================ Fig 10
def _sweep(a):
    kind, v, primed = a; prof = dict(BEST)
    if kind == 'alpha': prof['alpha'] = v
    elif kind == 'pT': prof['p_T'] = v * INP['p_W']
    elif kind == 'bT': prof['beta_T'] = v * INP['beta_W']
    elif kind == 'dose': prof['tip_dose'] = v
    return tstar(prof, {'primed': primed})
def fig10():
    grids = {'alpha': list(np.arange(0.02, 0.6001, 0.025)), 'pT': [1, 2, 3, 5, 7, 10], 'bT': list(np.logspace(math.log10(0.5), math.log10(8), 9)), 'dose': list(np.logspace(5, 9, 17))}   # grids: alpha 0.02-0.6 step 0.025, p_T ratios 1-10, beta_T ratios 0.5-8, dose 1e5-1e9
    meta = [('alpha', '(a) Hijacking factor', 'hijacking factor *α*', False), ('pT', '(b) Production ratio', 'production ratio *p*_{T}/*p*_{W}', False), ('bT', '(c) Entry ratio', 'entry ratio *β*_{T}/*β*_{W}', True), ('dose', '(d) TIP dose', 'TIP dose (copies/mL)', True)]
    F = Fig('fig10', 1290, 900); F.text(14, 30, 'Persistence threshold against TIP design parameters', mf.S_TITLE, bold=True, kind='title'); rows = []
    for k, (key, ttl, xl, lg) in enumerate(meta):
        X, Y = 8 + 640 * (k % 2), 46 + 400 * (k // 2); ser = []
        for primed, dash in ((1, None), (0, '8 5')):
            y = pmap(_sweep, [(key, v, primed) for v in grids[key]]); rows += [(key, primed, v, c) for v, c in zip(grids[key], y)]
            ser += segments(grids[key], y, BLUE, dash)
            if key == 'alpha': print('fig10 alpha primed' if primed else 'fig10 alpha abortive', [(round(v, 3), round(c, 3)) for v, c in zip(grids[key], y) if c == c][:5])
        xa = dict(label=xl, log=True) if lg else dict(label=xl)
        panel(F, X, Y, 630, 390, ttl, xa, dict(label='threshold *c*_{T}^{*} (d^{−1})', log=True), ser, bands=[(CW_BAND[0], CW_BAND[1], BAND_COL, '*c*_{W}')])
    mf.write_csv('fig10_threshold_vs_design', ['parameter', 'primed', 'value', 'cT_star'], rows)
    legend(F, 14, 860, [('primed cells tracked', BLUE, None), ('abortive-entry variant', BLUE, '8 5'), ('WT clearance range 5–15 d^{−1}', BAND_COL, None, 'patch')], 1200)
    return finish(F, 'fig10_threshold_vs_design.svg')

# ================================================================================================ Fig 11
SENS = [('*λ*', ['lambda']), ('*β*_{W}', ['beta_W']), ('*p*_{W}', ['p_W']), ('*δ*_{W} = *δ*_{Co}', ['delta_W', 'delta_Co']), ('*c*_{W}', ['c_W']), ('*d*_{T}', ['d_T']), ('*β*_{T}', ['beta_T']),
        ('TIP dose', ['tip_dose']), ('*p*_{T}', ['p_T']), ('*α*', ['alpha']), ('*k*_{W} = *k*_{Co}', ['k_W', 'k_Co']), ('*δ*_{N}', ['delta_N']), ('*K*_{F}', ['K_F']), ('*V*_{W}(0)', ['V_W0']), ('dosing time', ['tip_time'])]   # the fifteen parameters of Section 10.5
def _sens(a):
    names, f = a; prof = dict(BEST)
    for n in names: prof[n] = (BEST[n] if n in BEST else INP[n]) * f
    return tstar(prof)
def fig11():
    base = tstar(BEST); r = pmap(_sens, [(n, f) for _, n in SENS for f in (0.75, 1.25)])
    res = [(lab, (r[2 * i] / base - 1) * 100, (r[2 * i + 1] / base - 1) * 100) for i, (lab, _) in enumerate(SENS)]
    mf.write_csv('fig11_local_sensitivity', ['parameter', 'pct_change_cT_star_for_minus25', 'pct_change_for_plus25'], [(plain(l), a, b) for l, a, b in res])
    res.sort(key=lambda t: -max(abs(t[1]), abs(t[2]))); n = len(res); rh = 34
    W, H, L = 1000, 150 + n * rh, 250; pw = W - L - 40
    lim = mf.nice_step(0, max(max(abs(a), abs(b)) for _, a, b in res) * 1.1, 5) * 5; lim = math.ceil(max(max(abs(a), abs(b)) for _, a, b in res) * 1.1 / 20) * 20
    F = Fig('fig11', W, H); F.text(14, 30, f'Local sensitivity of the persistence threshold (baseline c_{{T}}^{{*}} = {base:.3f} d^{{−1}})', mf.S_TITLE, bold=True, kind='title')
    top = 60; xp = lambda v: L + (v + lim) / (2 * lim) * pw
    F.add(f"<rect x='{L}' y='{top}' width='{pw}' height='{n * rh}' fill='#fafafa'/>")
    for v in range(-lim, lim + 1, 20):
        F.line(xp(v), top, xp(v), top + n * rh, '#e4e4e4' if v else '#555', 1 if v else 1.6, check=False); F.line(xp(v), top + n * rh, xp(v), top + n * rh + 5, '#444', 1.2, check=False)
        F.text(xp(v), top + n * rh + 21, ('+' if v > 0 else '−' if v < 0 else '') + str(abs(v)), 13, 'middle', fill='#333', kind='tick')
    for i, (lab, a, b) in enumerate(res):
        y0 = top + i * rh; F.text(L - 14, y0 + rh / 2 + 5, lab, 14, 'end', fill='#111', kind='tick')
        for k, (val, col) in enumerate(((a, BLUE), (b, CORAL))):
            F.add(f"<rect x='{min(xp(0), xp(val)):.1f}' y='{y0 + 5 + k * 12:.1f}' width='{abs(xp(val) - xp(0)):.1f}' height='11' fill='{col}'/>")
    F.text(L + pw / 2, top + n * rh + 48, 'change in *c*_{T}^{*} (%)', mf.S_LABEL, 'middle', fill='#333', kind='axis')
    legend(F, L, top + n * rh + 76, [('parameter −25%', BLUE, None, 'patch'), ('parameter +25%', CORAL, None, 'patch')], 700)
    for lab, a, b in res[:5]: print(f'fig11 {plain(lab):14s} -25% -> {a:+.0f}%   +25% -> {b:+.0f}%')
    print(f'fig11 baseline c_T* = {base:.4f}; largest |c_T* after +-25%| = {max(max(a, b) for _, a, b in res) / 100 * base + base:.3f} (document: none above 0.4)')
    return finish(F, 'fig11_local_sensitivity.svg')
def plain(s): return rich(s, 14)[1]

# ================================================================================================ Fig 12
def _rho(r):
    ct = tstar(BEST, {'rho_T': r}, lo=0.01, hi=3.0); s = run({**BEST, 'rho_T': r, 'c_T': ct}); return ct, s['peakC']
def fig12():
    grid = list(np.linspace(0, 3e-3, 13)); r = pmap(_rho, grid)                                          # 13 values from 0 to 3e-3
    mf.write_csv('fig12_cytokine_tradeoff', ['rho_T', 'cT_star', 'peak_C_at_threshold'], [(g, a, b) for g, (a, b) in zip(grid, r)])
    F = Fig('fig12', 880, 560); F.text(14, 30, 'Trade-off between persistence threshold and peak cytokine', mf.S_TITLE, bold=True, kind='title')
    X, Y, W, H, RM = 8, 46, 860, 420, 150
    px, py = panel(F, X, Y, W, H, 'Innate stimulation of cytokine by free TIP', dict(label='TIP innate stimulation *ρ*_{T} (pg mL^{−1} d^{−1} per copy mL^{−1})', lo=0, hi=3e-3),
                   dict(label='threshold *c*_{T}^{*} (d^{−1})', lo=0, hi=0.5), [dict(x=grid, y=[a for a, _ in r], color=BLUE), dict(x=grid, y=[a for a, _ in r], color=BLUE, noline=True, marker='dot')], R=RM)
    L, T, B = 96, 36, 56; pw, ph = W - L - RM, H - T - B; pc = [b for _, b in r]
    lo, hi = math.floor(math.log10(min(pc))), math.ceil(math.log10(max(pc + [500])))
    yr = lambda v: Y + T + (1 - (math.log10(v) - lo) / (hi - lo)) * ph; xr = X + L + pw
    F.line(xr, Y + T, xr, Y + T + ph, RED, 1.6, check=False)
    for k in range(lo, hi + 1):
        yp = Y + T + (1 - (k - lo) / (hi - lo)) * ph; F.line(xr, yp, xr + 5, yp, RED, 1.2, check=False)
        m, pl, w = mf.pow10(k); F.add(f"<text x='{xr + 9:.1f}' y='{yp + 4.5:.1f}' font-size='13' text-anchor='start' fill='{RED}'>{m}</text>"); F.texts.append(dict(bb=(xr + 9, yp - 6, xr + 9 + w, yp + 8.5), plain=pl, kind='tick', note=False, owner=None, w=w))
    F.text(xr + 118, Y + T + ph / 2 + 5, 'peak cytokine at the threshold (pg/mL)', mf.S_LABEL, 'middle', fill=RED, kind='axis', rot=-90)
    pts = [(px(g), yr(c)) for g, c in zip(grid, pc)]
    F.add("<polyline fill='none' stroke='%s' stroke-width='2.4' stroke-linejoin='round' points='%s'/>" % (RED, ' '.join(f'{a:.1f},{b:.1f}' for a, b in pts)))
    for (a, b), (c, d) in zip(pts, pts[1:]): F.data.append((a, b, c, d))
    for a, b in pts: F.add(f"<circle cx='{a:.1f}' cy='{b:.1f}' r='3.8' fill='{RED}'/>")
    F.line(X + L, yr(500), xr, yr(500), RED, 1.5, '2 4', check=False); F.text(xr - 8, yr(500) + 22, '*C*_{crit} = 500 pg/mL', 13, 'end', fill=RED, kind='note', note=True)
    legend(F, 14, 500, [('threshold *c*_{T}^{*} (left axis)', BLUE, None), ('peak cytokine at the threshold (right axis)', RED, None), ('cytokine bound', RED, '2 4')], 840)
    print(f"fig12 c_T*: {r[0][0]:.3f} at rho_T=0 -> {r[-1][0]:.3f} at 3e-3 (document 0.256 -> 0.373); peak C {r[0][1]:.3g} -> {r[-1][1]:.3g} (document 5 -> 1.7e4)")
    return finish(F, 'fig12_cytokine_tradeoff.svg')

# ================================================================================================ Fig 14
def _vacc(a):
    A, xi, add = a; return tstar(BEST, {'A0': A, 'xi_T': xi, 'add_tip_clr': add}, adaptive=1)
def fig14():
    grid = [0.0] + [2.5 * i for i in range(1, 41)]
    cases = [('*ξ*_{T} = 0 (TIP invisible)', 0.0, 0, BLUE, None), ('*ξ*_{T} = 0.5', 0.5, 0, GREEN, None), ('*ξ*_{T} = 1, multiplicative', 1.0, 0, CORAL, None), ('*ξ*_{T} = 1, additive (worst case)', 1.0, 1, PURPLE, '8 5')]
    ser, leg, rows = [], [], []; naive = None
    for lab, xi, add, col, dash in cases:
        y = pmap(_vacc, [(A, xi, add) for A in grid]); rows += [(lab, A, c) for A, c in zip(grid, y)]
        ser += segments(grid, y, col, dash); ser.append(dict(x=[A for A, c in zip(grid, y) if c != c], y=[0.003] * sum(c != c for c in y), color=col, noline=True, marker='x')); leg.append((lab, col, dash))
        if xi == 1.0 and add == 0: naive = y[0]
        print(f"fig14 {plain(lab):34s} " + ', '.join(f'A={A:g}: {c:.3g}' for A, c in zip(grid, y) if A in (0, 12.5, 25, 50, 75) and c == c) + f"  | no extinction at {sum(c != c for c in y)} of {len(y)} titres")
    mf.write_csv('fig14a_cT_star_vs_antibody', ['assumption', 'A_vac', 'cT_star'], rows)
    Ag = [float(a) for a in np.linspace(0, 100, 101)]; reff = [run({'A0': a}, adaptive=1)['Reff0'] for a in Ag]
    def bis(f, lo=0.0, hi=100.0):
        for _ in range(30):
            m = (lo + hi) / 2
            if f(m): hi = m
            else: lo = m
        return hi
    Acrit = bis(lambda a: run({'A0': a}, adaptive=1)['Reff0'] <= 1.0)
    A60 = bis(lambda a: run({'A0': a, 'tip_dose': 0, 't_end': 60}, adaptive=1)['extinct'] == 1)
    mf.write_csv('fig14b_Reff_vs_antibody', ['A_vac', 'Reff'], list(zip(Ag, reff)))
    F = Fig('fig14', 1290, 600); F.text(14, 30, 'Vaccination and the persistence threshold', mf.S_TITLE, bold=True, kind='title')
    panel(F, 8, 46, 640, 440, '(a) Threshold against vaccine-induced antibody', dict(label='vaccine-induced antibody *A*_{vac} (µg/mL)', lo=0, hi=100), dict(label='threshold *c*_{T}^{*} (d^{−1})', log=True), ser,
          hlines=[(naive, 'naive', '#444', '2 4')], bands=[(CW_BAND[0], CW_BAND[1], BAND_COL, '*c*_{W}')], R=90)
    panel(F, 650, 46, 640, 440, '(b) Vaccine alone', dict(label='vaccine-induced antibody *A*_{vac} (µg/mL)', lo=0, hi=100), dict(label='*R*_{eff}(*A*)'), [dict(x=Ag, y=reff, color=BLUE)],
          hlines=[(1.0, '*R*_{eff} = 1', '#444', '2 4')], vlines=[(Acrit, f'*A* = {Acrit:.1f}', RED, '2 4'), (A60, f'*A* = {A60:.0f}', RED, '2 4')], R=100)
    legend(F, 14, 520, leg + [('no extinction for any *c*_{T} ≥ 0.003', '#444', None, 'x'), ('naive-host threshold', '#444', '2 4'), ('WT clearance range 5–15 d^{−1}', BAND_COL, None, 'patch')], 1250)
    print(f'fig14 A_crit (R_eff = 1) = {Acrit:.1f} (document 68.6); A for 60-day extinction without TIP = {A60:.1f} (document 94)')
    return finish(F, 'fig14_vaccination_threshold.svg')

# ================================================================================================ Fig 16
MCSETS = [(500, 2, 10, 8.73, 12000), (500, 2, 10, 20, 6000), (50, 2, 10, 8.73, 6000), (50, 2, 10, 4, 6000), (50, 2, 10, 2, 12000), (20, 2, 10, 2, 12000), (20, 2, 10, 1.5, 12000), (20, 2, 10, 1.2, 20000)]   # Table 16
def _mc(a):
    i, (p, dl, c, R, n) = a; AS.rng = np.random.default_rng(20261109 + i); return AS.mc(p, dl, c, R, n)
def fig16():
    mcs = pmap(_mc, list(enumerate(MCSETS)))
    th = [1 - AS.q_theory(p, dl, R) for p, dl, c, R, n in MCSETS]; se = [math.sqrt(f * (1 - f) / n) for f, (_, _, _, _, n) in zip(mcs, MCSETS)]
    mf.write_csv('fig16a_establishment', ['p', 'delta', 'c', 'Reff', 'lineages', 'theory_1_minus_q', 'monte_carlo', 'SE'], [(*s[:5], t, m, e) for s, t, m, e in zip(MCSETS, th, mcs, se)])
    F = Fig('fig16', 1290, 560); F.text(14, 30, 'Establishment of infection', mf.S_TITLE, bold=True, kind='title')
    X, Y, W, H = 8, 46, 640, 440
    px, py = panel(F, X, Y, W, H, '(a) Monte Carlo against theory', dict(label='theory 1 − *q* (Eq 26)', lo=0, hi=0.25), dict(label='Monte Carlo fraction that establish', lo=0, hi=0.25),
                   [dict(x=[0, 0.25], y=[0, 0.25], color=GREY, dash='8 5', width=1.8), dict(x=th, y=mcs, color=BLUE, noline=True, marker='dot')])
    for t, m, e in zip(th, mcs, se): F.line(px(t), py(m + 2 * e), px(t), py(m - 2 * e), BLUE, 1.6, check=False)
    ed = INP['delta_W'] + INP['delta_N'] * INP['C_rest'] / (INP['K_F'] + INP['C_rest'])                       # delta_eff = 2.05 d^-1 (Section 12.1)
    cases = [('naive host (*R*_{eff} = %.2f)' % run({})['R0'], run({})['R0'], BLUE, None), ('vaccinated host (*A* = 50 µg/mL, *R*_{eff} = %.2f)' % run({'A0': 50}, 1)['Reff0'], run({'A0': 50}, 1)['Reff0'], GREEN, None),
             ('*R*_{eff} = 3', 3.0, CORAL, None), ('*R*_{eff} = 0.9', 0.9, GREY, '8 5')]
    N = list(np.logspace(0, 3, 61)); ser = []; rows = []
    for lab, R, col, dash in cases:
        q = AS.q_theory(INP['p_W'], ed, R); pr = [1 - q ** n for n in N]; ser.append(dict(x=N, y=pr, color=col, dash=dash)); rows += [(plain(lab), n, v) for n, v in zip(N, pr)]
        print(f"fig16b {plain(lab):45s} q={q:.4f}  P(10)={1 - q ** 10:.3f}  P(100)={1 - q ** 100:.3f}  P(1000)={1 - q ** 1000:.3f}")
    mf.write_csv('fig16b_inoculum', ['case', 'inoculum', 'P_establish'], rows)
    panel(F, 650, 46, 640, 440, '(b) Inoculum size', dict(label='inoculum *N*_{0} (virions)', log=True), dict(label='probability of establishing infection', lo=0, hi=1), ser)
    legend(F, 14, 520, [('Monte Carlo ± 2 SE', BLUE, None, 'dot'), ('equality', GREY, '8 5')] + [(l, c, d) for l, _, c, d in cases], 1250)
    for s, t, m, e in zip(MCSETS, th, mcs, se): print(f'fig16a p={s[0]} R={s[3]:<5} theory {t:.4f}  MC {m:.4f} +- {e:.4f} ({s[4]} lineages)')
    return finish(F, 'fig16_establishment_probability.svg')

# ================================================================================================ Fig 17
TAIL_CT = [0.10, 0.20, 0.25, 0.26, 0.28, 0.30, 0.35, 0.40, 0.50, 0.70, 1.00]                      # the eleven values of Table 17
def _tail(a):
    i, ct = a; AS.rng = np.random.default_rng(20261109 + i)
    _, c = mf.rpt_traj({**BEST, 'c_T': ct, 'ext_floor': 0, 't_end': 30, 'dt_out': 0.01})
    env = np.array([c[k] for k in ('t', 'T', 'EW', 'IW', 'EC', 'IC', 'P', 'VW', 'VT', 'C')]).T
    pars = dict(bW=INP['beta_W'], bT=BEST['beta_T'], kW=INP['k_W'], kCo=INP['k_Co'], dW=INP['delta_W'], dCo=INP['delta_Co'], pW=INP['p_W'], alpha=BEST['alpha'], cW=INP['c_W'], dN=INP['delta_N'], KF=INP['K_F'])
    ext, esc, t0, tot = AS.tail_sim(env, pars, 300, 200); det = run({**BEST, 'c_T': ct})['extinct']
    return ct, ext, esc, t0, tot, det
def wilson(k, n, z=1.96):
    p = k / n; d = 1 + z * z / n; c = (p + z * z / (2 * n)) / d; h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d; return c - h, c + h
def fig17():
    r = pmap(_tail, list(enumerate(TAIL_CT))); ts = tstar(BEST)
    mf.write_csv('fig17_stochastic_tail', ['c_T', 'extinct_of_200', 'escaped', 'tail_start_day', 'tail_start_copies', 'deterministic_extinct'], r)
    F = Fig('fig17', 880, 560); F.text(14, 30, 'Extinction probability in the stochastic tail', mf.S_TITLE, bold=True, kind='title')
    x = [a[0] for a in r]; p = [a[1] / 200 for a in r]; ci = [wilson(a[1], 200) for a in r]
    det = dict(x=[min(x), ts, ts, max(x)], y=[1, 1, 0, 0], color=GREY, dash='8 5', width=1.8)
    px, py = panel(F, 8, 46, 860, 420, '(a) 200 replicates per point', dict(label=CT_LAB, log=True, lo=math.log10(0.1), hi=0), dict(label='probability of WT extinction', lo=0, hi=1.0),
                   [det, dict(x=x, y=p, color=BLUE), dict(x=x, y=p, color=BLUE, noline=True, marker='dot')], vlines=[(ts, '', RED, '6 4')])
    for xv, (lo, hi) in zip(x, ci): F.line(px(xv), py(hi), px(xv), py(lo), BLUE, 1.6, check=False)
    legend(F, 14, 500, [('stochastic tail (95% Wilson interval)', BLUE, None, 'dot'), ('deterministic outcome', GREY, '8 5'), (f'deterministic threshold *c*_{{T}}^{{*}} = {ts:.3f} d^{{−1}}', RED, '6 4')], 840)
    for a in r: print(f'fig17 c_T={a[0]:<5} extinct {a[1]:>3}/200 (P={a[1] / 200:.3f}) escaped {a[2]:>3} tail start day {a[3]:.1f}, {a[4]} copies; deterministic extinct={int(a[5])}')
    print(f'fig17 deterministic threshold c_T* = {ts:.4f} d^-1')
    return finish(F, 'fig17_stochastic_tail.svg')

ALL = dict(fig2=fig2, fig4=fig4, fig6=fig6, fig8=fig8, fig9=fig9, fig10=fig10, fig11=fig11, fig12=fig12, fig14=fig14, fig16=fig16, fig17=fig17)
if __name__ == '__main__':
    want = [a for a in sys.argv[2:] if a in ALL] or list(ALL)
    bad = sum(ALL[k]() for k in want)
    print('make_figures_analysis: all layout checks passed' if bad == 0 else f'make_figures_analysis: {bad} layout problem(s)')
    sys.exit(1 if bad else 0)
