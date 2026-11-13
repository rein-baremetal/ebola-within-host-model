#!/usr/bin/env python3
# verify_figures.py - confirms that the numbers drawn in the engine-written SVGs equal the engine's own output.
# It reads the SVG geometry (tick marks + tick labels give the axis calibration), converts every plotted point back to data
# units and compares it with the engine's text/CSV output.   Usage:
#   ./engine input.txt sample_trajectory.txt errors.txt         (writes sample_trajectory.svg)
#   python3 verify_figures.py sample_trajectory.svg sample_trajectory.txt [untreated_output.txt]
#   python3 verify_figures.py --pareto screen_pareto.svg screen_samples.csv
import re, sys, math, csv
if len(sys.argv) < 3:
    sys.exit(__doc__ if __doc__ else 'usage: verify_figures.py sample_trajectory.svg sample_trajectory.txt untreated_output.txt | --pareto screen_pareto.svg screen_samples.csv')
def num(label):
    m = re.match(r"10<tspan[^>]*>([\u2212-]?\d+)", label)
    if m: return ('log', int(m.group(1).replace('\u2212', '-')))
    label = re.sub(r'<[^>]*>', '', label)
    return ('lin', float(label.replace(',', '').replace('\u2212', '-')))
def calib(pairs):            # pairs (pixel, label) -> function pixel -> data value
    (p0, l0), (p1, l1) = pairs[0], pairs[-1]; k0, v0 = num(l0); k1, v1 = num(l1)
    f = lambda p: v0 + (p - p0) / (p1 - p0) * (v1 - v0)
    return (lambda p: 10 ** f(p)) if k0 == 'log' else f
def axes(chunk):
    yl = [(float(m.group(2)), m.group(3)) for m in re.finditer(r"<text x='([\d.]+)' y='([\d.]+)' font-size='13' text-anchor='end' fill='#333'>(.*?)</text>", chunk)]
    xl = [(float(m.group(1)), m.group(3)) for m in re.finditer(r"<text x='([\d.]+)' y='([\d.]+)' font-size='13' text-anchor='middle' fill='#333'>(.*?)</text>", chunk)]
    yl = [(y - 4.5, l) for y, l in yl]
    return xl, yl
def traj(svg, out, out_u):
    t = open(svg, encoding='utf8').read()
    chunks = re.split(r"(?=<text [^>]*font-weight='bold'[^>]*>(?:<tspan[^>]*>)?\([abcd]\))", t)[1:]
    def cols(fn):
        rows = [l.split() for l in open(fn) if l.strip() and not l.startswith('#')]
        return {k: [float(r[i]) for r in rows] for i, k in enumerate('t T EW IW VW EC IC VT C H J'.split())}
    A = cols(out); U = cols(out_u) if out_u else None
    spec = [[('VW', A), ('VT', A)] + ([('VW', U)] if U else []), [('IW', A), ('IC', A)] + ([('IW', U)] if U else []), [('C', A)] + ([('C', U)] if U else []), [('H', A)] + ([('H', U)] if U else [])]
    worst = 0.0; n = 0
    for ch, sp in zip(chunks, spec):
        xl, yl = axes(ch); fx, fy = calib(xl), calib(yl)
        pls = re.findall(r"<polyline[^>]*points='([^']*)'", ch)
        assert len(pls) == len(sp), (len(pls), len(sp))
        logy = num(yl[0][1])[0] == 'log'; rng = abs(fy(yl[0][0]) - fy(yl[-1][0])) if not logy else 1.0
        for pl, (key, D) in zip(pls, sp):
            pts = [tuple(map(float, p.split(','))) for p in pl.split()]
            assert len(pts) == len(D['t'])
            for (px, py), tt, v in zip(pts, D['t'], D[key]):
                got = fy(py); exp = v
                if logy: err = abs(math.log10(got) - math.log10(max(exp, 1e-4)))
                else: err = abs(got - exp) / rng
                worst = max(worst, err, abs(fx(px) - tt) / 30); n += 1
    print(f'{svg}: {n} plotted points compared with {out}: worst deviation {worst:.2e} (fraction of axis range; log axes in decades)')
    return worst < 2e-3
def pareto(svg, csvf):
    t = open(svg, encoding='utf8').read(); xl, yl = axes(t); fx, fy = calib(xl), calib(yl)
    cs = [(float(m.group(1)), float(m.group(2))) for m in re.finditer(r"<circle cx='([\d.]+)' cy='([\d.]+)' r='(?:3|2\.4|4\.8)'", t)]
    ybot = max(p for p, _ in yl); cs = [c for c in cs if c[1] <= ybot + 1]      # drop the legend markers below the plot
    got = sorted((round(math.log10(fx(a)), 3), round(fy(b), 1)) for a, b in cs)
    rows = list(csv.DictReader(open(csvf)))
    exp = sorted((round(math.log10(max(float(r['peakV']), 1e-4)), 3), round(min(float(r['peakC']), fy(yl[-1][0])), 1)) for r in rows if True)
    # points above the axis top would be clamped; the engine plots every valid sample
    import bisect
    ex = sorted(exp); worst = 0.0
    for g in got:                                   # nearest engine sample for every drawn point (log10 peak V, peak C)
        k = bisect.bisect_left(ex, g); cand = ex[max(0, k - 40):k + 40]
        worst = max(worst, min(max(abs(g[0] - e[0]) / 0.02, abs(g[1] - e[1]) / 1.0) for e in cand))
    ok = len(got) == len(exp) and worst <= 1.0
    print(f'{svg}: {len(got)} circles vs {len(exp)} samples in {csvf}: {"match" if ok else "MISMATCH"}'); return ok
if __name__ == '__main__':
    a = sys.argv[1:]
    ok = pareto(a[1], a[2]) if a and a[0] == '--pareto' else traj(a[0], a[1], a[2] if len(a) > 2 else None)
    sys.exit(0 if ok else 1)
