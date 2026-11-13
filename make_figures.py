#!/usr/bin/env python3
# make_figures.py - figures that code_annotated.cpp does not draw itself (Python 3 + numpy; the C++ engine stays dependency-free).
#
# Every plotted number comes from running the package's own programs:
#   * ./rpt_adapt  (adaptive-immunity engine, same equations as code_annotated.cpp; built from rpt_adapt.cpp)
#   * analysis_scripts.py : run(), BEST (screened-best TIP profile)  -> used for Fig 13 and Fig 15
# Nothing is typed in by hand: each figure also writes the exact numbers it plots to figure_data/<name>.csv.
#
#   Fig 5  fig05_numerical_verification.svg   (a) RK4 convergence, (b) return to the infection-free equilibrium
#   Fig 7  fig07_model_schematic.svg          structure of the extended model (no numbers; drawn from the right-hand sides in rpt_adapt.cpp)
#   Fig 13 fig13_vaccination_trajectories.svg naive / vaccinated hosts, untreated / with TIP
#   Fig 15 fig15_vaccine_antibody_sweep.svg   peak virus and peak cytokine against vaccine-induced antibody
#
# Usage:  g++ -std=c++20 -O2 rpt_adapt.cpp -o rpt_adapt        (once)
#         python3 make_figures.py [outdir]                       (outdir defaults to the package folder)
# Style rules, text-overlap / box-overflow checks: see the Fig class below; any violation is printed and the exit code is 1.
import csv, importlib.util, math, os, subprocess, sys, tempfile
from PIL import ImageFont

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = sys.argv[1] if len(sys.argv) > 1 else HERE
os.makedirs(os.path.join(OUT, 'figure_data'), exist_ok=True)

# ---- reuse the package's own driver functions (run, BEST) -------------------------------------------
_spec = importlib.util.spec_from_file_location('analysis_scripts', os.path.join(HERE, 'analysis_scripts.py'))
AS = importlib.util.module_from_spec(_spec); _spec.loader.exec_module(AS)
EXE = AS.EXE
if not os.path.exists(EXE):
    sys.exit('rpt_adapt not built: g++ -std=c++20 -O2 rpt_adapt.cpp -o rpt_adapt')

# ================================================================================================
# Text measurement and mini-markup:  *V*  italic,  _{W} subscript,  ^{2} superscript   (same notation as the C++ engine)
# ================================================================================================
_FD = '/usr/share/fonts/truetype/dejavu/'
_FONTS = {}
def _font(size, bold=False, ital=False):
    key = (round(size * 4), bold, ital)
    if key not in _FONTS:
        name = 'DejaVuSans' + ('-Bold' if bold else '') + ('-Oblique' if ital and not bold else '') + '.ttf'
        try: _FONTS[key] = ImageFont.truetype(_FD + name, size)
        except OSError: _FONTS[key] = None
    return _FONTS[key]
def text_width(s, size, bold=False, ital=False):
    f = _font(size, bold, ital)
    return f.getlength(s) if f else 0.6 * size * len(s)      # DejaVu Sans is wider than Arial: the estimate is conservative

def rich(s, fs, bold=False):
    """-> (svg tspans, plain text, width px).  Sub/superscripts are real (smaller, shifted) and shifts are undone afterwards."""
    segs, cur, ital, shift, i = [], '', False, 0.0, 0
    def flush(sz=100, sh=0.0, it=False):
        nonlocal cur
        if cur: segs.append((cur, sz, sh, it)); cur = ''
    while i < len(s):
        c = s[i]
        if c == '*': flush(100, 0, ital); ital = not ital; i += 1; continue
        if c in '_^' and i + 1 < len(s) and s[i + 1] == '{':
            flush(100, 0, ital); e = s.index('}', i)
            segs.append((s[i + 2:e], 75, (0.22 if c == '_' else -0.38) * fs, ital)); i = e + 1; continue
        cur += c; i += 1
    flush(100, 0, ital)
    out, plain, w, prev = '', '', 0.0, 0.0
    for t, sz, sh, it in segs:
        esc = t.replace('&', '&amp;').replace('<', '&lt;')
        a = ''
        if it: a += " font-style='italic'"
        if sz < 100: a += f" font-size='{sz}%'"
        if abs(sh - prev) > 1e-9: a += f" dy='{sh - prev:.2f}'"
        out += f'<tspan{a}>{esc}</tspan>'; plain += t; prev = sh
        w += text_width(t, fs * sz / 100, bold, it)
    if prev: out += f"<tspan dy='{-prev:.2f}'></tspan>"
    return out, plain, w

# ================================================================================================
# Geometry helpers for the checks
# ================================================================================================
def _seg_seg_dist(a, b):
    (x1, y1, x2, y2), (x3, y3, x4, y4) = a, b
    def cross(ax, ay, bx, by): return ax * by - ay * bx
    d1 = cross(x2 - x1, y2 - y1, x3 - x1, y3 - y1); d2 = cross(x2 - x1, y2 - y1, x4 - x1, y4 - y1)
    d3 = cross(x4 - x3, y4 - y3, x1 - x3, y1 - y3); d4 = cross(x4 - x3, y4 - y3, x2 - x3, y2 - y3)
    if d1 * d2 < 0 and d3 * d4 < 0: return 0.0
    def pd(px, py, sx1, sy1, sx2, sy2):
        dx, dy = sx2 - sx1, sy2 - sy1; L2 = dx * dx + dy * dy
        t = 0 if L2 == 0 else max(0, min(1, ((px - sx1) * dx + (py - sy1) * dy) / L2))
        return math.hypot(px - (sx1 + t * dx), py - (sy1 + t * dy))
    return min(pd(x1, y1, *b), pd(x2, y2, *b), pd(x3, y3, *a), pd(x4, y4, *a))
def seg_rect_dist(seg, r):
    x0, y0, x1, y1 = r
    inside = lambda px, py: x0 < px < x1 and y0 < py < y1
    if inside(seg[0], seg[1]) or inside(seg[2], seg[3]): return 0.0
    return min(_seg_seg_dist(seg, e) for e in ((x0, y0, x1, y0), (x1, y0, x1, y1), (x1, y1, x0, y1), (x0, y1, x0, y0)))
def rect_rect_dist(a, b):
    dx = max(b[0] - a[2], a[0] - b[2], 0); dy = max(b[1] - a[3], a[1] - b[3], 0); return math.hypot(dx, dy)

# ================================================================================================
# Fig: one SVG canvas that records every text box, box, arrow and data line so the layout can be checked.
# ================================================================================================
FONT = 'Arial,Helvetica,sans-serif'
S_TITLE, S_PTITLE, S_BOX, S_LABEL, S_TICK, S_LEG, S_NOTE = 21, 17, 19, 16, 13, 14, 13
class Fig:
    def __init__(self, name, W, H):
        self.name, self.W, self.H = name, W, H
        self.parts, self.texts, self.boxes, self.arrows, self.lines, self.data, self.problems = [], [], [], [], [], [], []
    def add(self, s): self.parts.append(s)
    # -- text: registers its bounding box (kind: 'label' = must keep 8 px from boxes/lines, 'tick', 'title', 'box' = text inside a box)
    def text(self, x, y, markup, size=S_LABEL, anchor='start', bold=False, fill='#111', kind='label', note=False, rot=None, owner=None):
        svg, plain, w = rich(markup, size, bold)
        x0 = x - w if anchor == 'end' else x - w / 2 if anchor == 'middle' else x
        a = f"x='{x:.1f}' y='{y:.1f}' font-size='{size}' text-anchor='{anchor}'" + (" font-weight='bold'" if bold else '') + f" fill='{fill}'"
        if rot is not None:
            self.add(f"<text {a} transform='rotate({rot} {x:.1f} {y:.1f})'>{svg}</text>")
            bb = (x - size, y - w / 2 - 2, x + 0.3 * size, y + w / 2 + 2) if rot == -90 else None
            # rotated label: x,y is the centre of the text line; box width is the text height
            bb = (x - 0.85 * size, y - w / 2, x + 0.3 * size, y + w / 2)
        else:
            self.add(f"<text {a}>{svg}</text>")
            bb = (x0, y - 0.80 * size, x0 + w, y + 0.25 * size)
        self.texts.append(dict(bb=bb, plain=plain, kind=kind, note=note, owner=owner, w=w))
        return bb
    def line(self, x0, y0, x1, y1, color='#222', sw=1.6, dash=None, check=True):
        self.add(f"<line x1='{x0:.1f}' y1='{y0:.1f}' x2='{x1:.1f}' y2='{y1:.1f}' stroke='{color}' stroke-width='{sw}'" + (f" stroke-dasharray='{dash}'" if dash else '') + '/>')
        if check: self.lines.append((x0, y0, x1, y1))
    def box(self, x, y, w, h, l1, l2=None, name=None):
        self.add(f"<rect x='{x}' y='{y}' width='{w}' height='{h}' fill='white' stroke='#000' stroke-width='1.6'/>")
        self.boxes.append(dict(r=(x, y, x + w, y + h), name=name or l1))
        cy = y + h / 2
        self.text(x + w / 2, cy - 2 if l2 else cy + 6, l1, S_BOX, 'middle', kind='box', owner=name or l1)
        if l2: self.text(x + w / 2, cy + 22, l2, S_LABEL, 'middle', kind='box', owner=name or l1, fill='#333')
    def frame(self, x0, y0, x1, y1):
        self.add(f"<rect x='{x0}' y='{y0}' width='{x1 - x0}' height='{y1 - y0}' fill='none' stroke='#9a9a9a' stroke-width='1.2'/>")
        for e in ((x0, y0, x1, y0), (x1, y0, x1, y1), (x1, y1, x0, y1), (x0, y1, x0, y0)): self.lines.append(e)
    def arrow(self, x0, y0, x1, y1, head='arrow', dash=False, name=None, joins=None, tail=False):
        """Straight horizontal/vertical connector. head: 'arrow' (filled head) | 'flat' (flat bar = inhibition) | None."""
        assert x0 == x1 or y0 == y1, 'arrows must be horizontal or vertical'
        ux, uy = (x1 - x0), (y1 - y0); L = math.hypot(ux, uy); ux, uy = ux / L, uy / L
        ex, ey = x1, y1
        if head == 'arrow': ex, ey = x1 - 11 * ux, y1 - 11 * uy
        self.line(x0, y0, ex, ey, '#000', 1.6, '7 5' if dash else None, check=False)
        if head == 'arrow':
            px, py = -uy, ux
            self.add(f"<polygon points='{x1:.1f},{y1:.1f} {ex + 5.5 * px:.1f},{ey + 5.5 * py:.1f} {ex - 5.5 * px:.1f},{ey - 5.5 * py:.1f}' fill='#000'/>")
        elif head == 'flat':
            px, py = -uy, ux
            self.add(f"<line x1='{x1 + 8 * px:.1f}' y1='{y1 + 8 * py:.1f}' x2='{x1 - 8 * px:.1f}' y2='{y1 - 8 * py:.1f}' stroke='#000' stroke-width='2.4'/>")
        self.arrows.append(dict(s=(x0, y0, x1, y1), name=name, joins=joins))
    # -- the checks ---------------------------------------------------------------------------
    def check(self):
        P = self.problems
        for i, t in enumerate(self.texts):
            x0, y0, x1, y1 = t['bb']
            if x0 < 0 or y0 < 0 or x1 > self.W or y1 > self.H: P.append(f"text '{t['plain']}' leaves the canvas")
            for u in self.texts[i + 1:]:
                if rect_rect_dist(t['bb'], u['bb']) == 0 and (min(t['bb'][2], u['bb'][2]) - max(t['bb'][0], u['bb'][0]) > 0.5) and (min(t['bb'][3], u['bb'][3]) - max(t['bb'][1], u['bb'][1]) > 0.5):
                    P.append(f"text '{t['plain']}' overlaps text '{u['plain']}'")
        for t in self.texts:                       # text overflow inside boxes: longest line + 20 px padding on both sides must fit
            if t['kind'] == 'box':
                b = next(b for b in self.boxes if b['name'] == t['owner'])
                if t['w'] + 40 > b['r'][2] - b['r'][0] + 1e-6: P.append(f"text '{t['plain']}' + 20 px padding overflows box '{b['name']}'")
        for t in self.texts:                       # labels: >= 8 px from every box, connector and frame line; box text must not touch connectors
            if t['kind'] == 'label':
                for b in self.boxes:
                    if rect_rect_dist(t['bb'], b['r']) < 8: P.append(f"label '{t['plain']}' is closer than 8 px to box '{b['name']}'")
                for a in self.arrows:
                    if seg_rect_dist(a['s'], t['bb']) < 8: P.append(f"label '{t['plain']}' is closer than 8 px to a connector ({a['name']})")
                for l in self.lines:
                    if seg_rect_dist(l, t['bb']) < 8: P.append(f"label '{t['plain']}' is closer than 8 px to a frame line")
            elif t['kind'] == 'box':
                for a in self.arrows:
                    if seg_rect_dist(a['s'], t['bb']) == 0: P.append(f"connector {a['name']} crosses box text '{t['plain']}'")
            if t['note']:                          # annotations on graphs must stay off the data lines (4 px clearance)
                for s in self.data:
                    if seg_rect_dist(s, t['bb']) < 4: P.append(f"annotation '{t['plain']}' lies on a data line"); break
        for i, a in enumerate(self.arrows):        # connectors never cross each other (a branch that ends ON another connector is declared with joins=)
            for b in self.arrows[i + 1:]:
                if b['name'] in (a['joins'] or ()) or a['name'] in (b['joins'] or ()): continue
                if _seg_seg_dist(a['s'], b['s']) < 1.0: P.append(f"connectors {a['name']} and {b['name']} touch or cross")
        for a in self.arrows:                      # connectors never pass through boxes (only touch the edge they start/end at)
            for b in self.boxes:
                x0, y0, x1, y1 = b['r']; r = (x0 + 1.5, y0 + 1.5, x1 - 1.5, y1 - 1.5)
                if seg_rect_dist(a['s'], r) == 0: P.append(f"connector {a['name']} enters box '{b['name']}'")
        print(f"layout check [{self.name}]: {len(P)} problem(s) in {len(self.texts)} text items, {len(self.arrows)} connectors, {len(self.boxes)} boxes")
        for p in P: print('  LAYOUT:', p)
        return len(P)
    def save(self, path):
        with open(path, 'w', encoding='utf8') as f:
            f.write(f"<svg xmlns='http://www.w3.org/2000/svg' width='{self.W}' height='{self.H}' viewBox='0 0 {self.W} {self.H}' font-family='{FONT}' xml:space='preserve'>\n")
            f.write("<rect width='100%' height='100%' fill='white'/>\n" + '\n'.join(self.parts) + '\n</svg>\n')

# ================================================================================================
# Graph panels (same look as panel() in code_annotated.cpp): bottom+left spines only, light grid, round-number ticks, no minor ticks
# ================================================================================================
def nice_step(lo, hi, target):
    raw = max(hi - lo, 1e-300) / target; mag = 10 ** math.floor(math.log10(raw)); r = raw / mag
    return (1 if r < 1.5 else 2 if r < 3 else 5 if r < 7 else 10) * mag
def fmt_tick(v, step):
    dec = max(0, math.ceil(-math.log10(step) - 1e-9)); s = f'{0.0 if abs(v) < step * 1e-9 else v:.{dec}f}'
    sign = '−' if s.startswith('-') else ''; s = s.lstrip('-')
    ip, _, fp = s.partition('.'); ip = f'{int(ip):,}'
    return sign + ip + ('.' + fp if fp else '')
def pow10(k, fs=S_TICK):
    e = ('−' if k < 0 else '') + str(abs(k))
    return f"10<tspan dy='{-0.38 * fs:.2f}' font-size='75%'>{e}</tspan><tspan dy='{0.38 * fs:.2f}'></tspan>", '10' + e, text_width('10', fs) + text_width(e, fs * 0.75)

def panel(F, X, Y, W, H, title, xax, yax, series, hlines=(), notes=(), clamp=1e-4, R=70, vlines=(), bands=()):
    """xax/yax: dict(label=, log=False, lo=None, hi=None).  series: dict(x=[], y=[], color=, dash=None, width=2.4).
    hlines: (value, markup label, colour, dash).  Their labels sit in the right margin, never on the data.
    notes: (x_data, y_data, markup, anchor) annotations placed by the caller in empty plot area (checked against the data)."""
    L, T, B = 96, 36, 56; pw, ph = W - L - R, H - T - B
    def tf(v, ax): return math.log10(max(v, clamp)) if ax.get('log') else v
    allx = [tf(v, xax) for s in series for v in s['x']]; ally = [tf(v, yax) for s in series for v in s['y']] + [tf(h[0], yax) for h in hlines] + [tf(v, yax) for b in bands for v in b[:2]]
    def axis(ax, vals, pad):
        lo, hi = (ax['lo'] if ax.get('lo') is not None else min(vals)), (ax['hi'] if ax.get('hi') is not None else max(vals))
        if ax.get('log'):
            lo, hi = math.floor(lo + 1e-9), math.ceil(hi - 1e-9)
            if hi - lo < 1: hi = lo + 1
            step = 2 if hi - lo > 6 else 1
        else:
            if ax.get('lo') is None: lo = min(0.0, lo)
            if ax.get('hi') is None: hi = hi * pad + 1e-12
            step = nice_step(lo, hi, 5)
            if ax.get('hi') is None: hi = math.ceil(hi / step - 1e-9) * step
            if ax.get('lo') is None: lo = math.floor(lo / step + 1e-9) * step
        return lo, hi, step
    x0, x1, xs = axis(xax, allx, 1.0); y0, y1, ys = axis(yax, ally, 1.06)
    px = lambda v: X + L + (tf(v, xax) - x0) / (x1 - x0) * pw
    py = lambda v: Y + T + (1 - (tf(v, yax) - y0) / (y1 - y0)) * ph
    F.text(X + L, Y + 22, title, S_PTITLE, bold=True, fill='#1a1a1a', kind='title')
    F.add(f"<rect x='{X + L}' y='{Y + T}' width='{pw}' height='{ph}' fill='#fafafa'/>")
    for lo_, hi_, col_, lab_ in bands:             # shaded reference band (label in the right margin)
        yt, yb = sorted((py(hi_), py(lo_))); yt, yb = max(Y + T, yt), min(Y + T + ph, yb)
        F.add(f"<rect x='{X + L}' y='{yt:.1f}' width='{pw}' height='{yb - yt:.1f}' fill='{col_}' fill-opacity='0.35'/>")
        F.text(X + L + pw + 8, (yt + yb) / 2 + 4.5, lab_, S_NOTE, fill='#555', kind='ref', note=True)
    def ticks(lo, hi, st):
        v = math.ceil(lo / st - 1e-9) * st
        while v <= hi + st * 1e-9: yield v; v += st
    for v in ticks(x0, x1, xs):
        xp = X + L + (v - x0) / (x1 - x0) * pw
        F.line(xp, Y + T, xp, Y + T + ph, '#e4e4e4', 1, check=False); F.line(xp, Y + T + ph, xp, Y + T + ph + 5, '#444', 1.2, check=False)
        if xax.get('log'):
            m, pl, w = pow10(int(round(v))); F.add(f"<text x='{xp:.1f}' y='{Y + T + ph + 21:.1f}' font-size='{S_TICK}' text-anchor='middle' fill='#333'>{m}</text>")
            F.texts.append(dict(bb=(xp - w / 2, Y + T + ph + 21 - 0.8 * S_TICK, xp + w / 2, Y + T + ph + 21 + 4), plain=pl, kind='tick', note=False, owner=None, w=w))
        else: F.text(xp, Y + T + ph + 21, fmt_tick(v, xs), S_TICK, 'middle', fill='#333', kind='tick')
    for v in ticks(y0, y1, ys):
        yp = Y + T + (1 - (v - y0) / (y1 - y0)) * ph
        F.line(X + L, yp, X + L + pw, yp, '#e4e4e4', 1, check=False); F.line(X + L - 5, yp, X + L, yp, '#444', 1.2, check=False)
        if yax.get('log'):
            m, pl, w = pow10(int(round(v))); F.add(f"<text x='{X + L - 9:.1f}' y='{yp + 4.5:.1f}' font-size='{S_TICK}' text-anchor='end' fill='#333'>{m}</text>")
            F.texts.append(dict(bb=(X + L - 9 - w, yp + 4.5 - 0.8 * S_TICK, X + L - 9, yp + 4.5 + 4), plain=pl, kind='tick', note=False, owner=None, w=w))
        else: F.text(X + L - 9, yp + 4.5, fmt_tick(v, ys), S_TICK, 'end', fill='#333', kind='tick')
    F.line(X + L, Y + T + ph, X + L + pw, Y + T + ph, '#222', 1.6, check=False); F.line(X + L, Y + T, X + L, Y + T + ph, '#222', 1.6, check=False)
    F.text(X + L + pw / 2, Y + H - 12, xax['label'], S_LABEL, 'middle', fill='#333', kind='axis')
    F.text(X + 18, Y + T + ph / 2 + 5, yax['label'], S_LABEL, 'middle', fill='#333', kind='axis', rot=-90)
    for val, lab, col, dash in hlines:             # reference lines; label in the right margin
        F.line(X + L, py(val), X + L + pw, py(val), col, 1.5, dash, check=False)
        F.text(X + L + pw + 8, py(val) + 4.5, lab, S_NOTE, fill=col, kind='ref', note=True)
    for v_, lab_, col_, dash_ in vlines:           # vertical reference lines; label just right of the line at the top of the plot
        F.line(px(v_), Y + T, px(v_), Y + T + ph, col_, 1.5, dash_, check=False)
        if lab_: F.text(px(v_) + 6, Y + T + 16, lab_, S_NOTE, fill=col_, kind='note', note=True)
    for s in series:
        pts = ' '.join(f'{px(a):.1f},{py(b):.1f}' for a, b in zip(s['x'], s['y']))
        if not s.get('noline'): F.add(f"<polyline fill='none' stroke='{s['color']}' stroke-width='{s.get('width', 2.4)}' stroke-linejoin='round'" + (f" stroke-dasharray='{s['dash']}'" if s.get('dash') else '') + f" points='{pts}'/>")
        for k in range(1, 1 if s.get('noline') else len(s['x'])): F.data.append((px(s['x'][k - 1]), py(s['y'][k - 1]), px(s['x'][k]), py(s['y'][k])))
    for s_ in series:                              # point markers (x = no result, tri = threshold, default dot)
        if s_.get('marker'):
            for a_, b_ in zip(s_['x'], s_['y']):
                cx_, cy_, c_ = px(a_), py(b_), s_['color']
                if s_['marker'] == 'x': F.add(f"<path d='M{cx_ - 5:.1f} {cy_ - 5:.1f}L{cx_ + 5:.1f} {cy_ + 5:.1f}M{cx_ - 5:.1f} {cy_ + 5:.1f}L{cx_ + 5:.1f} {cy_ - 5:.1f}' stroke='{c_}' stroke-width='2'/>")
                elif s_['marker'] == 'tri': F.add(f"<polygon points='{cx_:.1f},{cy_ - 6:.1f} {cx_ - 6:.1f},{cy_ + 5:.1f} {cx_ + 6:.1f},{cy_ + 5:.1f}' fill='{c_}'/>")
                else: F.add(f"<circle cx='{cx_:.1f}' cy='{cy_:.1f}' r='3.8' fill='{c_}'/>")
                F.data.append((cx_ - 5, cy_ - 5, cx_ + 5, cy_ + 5))
    for xv, yv, mk, anc in notes:
        F.text(px(xv), py(yv), mk, S_NOTE, anc, fill='#333', kind='note', note=True)
    return px, py

def legend(F, x, y, items, maxw):
    """items: (markup, colour, dash[, marker]).  One or more rows under the panels (an empty area)."""
    cx, cy = x, y
    for it in items:
        mk, col, dash = it[:3]; marker = it[3] if len(it) > 3 else None
        _, plain, w = rich(mk, S_LEG); iw = 30 + w
        if cx > x and cx + iw > x + maxw: cx, cy = x, cy + 26
        if marker == 'x': F.add(f"<path d='M{cx + 7:.1f} {cy - 5:.1f}L{cx + 17:.1f} {cy + 5:.1f}M{cx + 7:.1f} {cy + 5:.1f}L{cx + 17:.1f} {cy - 5:.1f}' stroke='{col}' stroke-width='2'/>")
        elif marker == 'tri': F.add(f"<polygon points='{cx + 12:.1f},{cy - 6:.1f} {cx + 6:.1f},{cy + 5:.1f} {cx + 18:.1f},{cy + 5:.1f}' fill='{col}'/>")
        elif marker == 'patch': F.add(f"<rect x='{cx}' y='{cy - 7}' width='24' height='14' fill='{col}' fill-opacity='0.55' stroke='#777' stroke-width='0.8'/>")
        else: F.line(cx, cy, cx + 24, cy, col, 2.4, dash, check=False)
        F.text(cx + 30, cy + 5, mk, S_LEG, fill='#222', kind='legend')
        cx += iw + 26
    return cy + 26

# ================================================================================================
# Data access: parameters from input.txt, runs through the package's own run() / rpt_adapt
# ================================================================================================
def read_input(path):
    d = {}
    for ln in open(path, encoding='utf8'):
        ln = ln.split('#')[0].strip()
        if '=' in ln:
            k, v = [t.strip() for t in ln.split('=', 1)]
            try: d[k] = float(v)
            except ValueError: pass
    return d
def rpt_traj(over):
    """Run rpt_adapt with the given overrides and return (summary dict, trajectory dict of column lists)."""
    with tempfile.TemporaryDirectory() as td:
        f = os.path.join(td, 'traj.csv'); s = AS.run(over, traj=f, adaptive=over.get('adaptive', 0)) if 'adaptive' not in over else AS.run(over, traj=f)
        rows = list(csv.reader(open(f))); hdr = rows[0]
        cols = {h: [float(r[i]) for r in rows[1:]] for i, h in enumerate(hdr)}
    return s, cols
def write_csv(name, header, rows):
    with open(os.path.join(OUT, 'figure_data', name + '.csv'), 'w', newline='') as f:
        w = csv.writer(f); w.writerow(header); w.writerows(rows)

# ================================================================================================
# Figure 5: numerical verification  (a) RK4 convergence  (b) return to the infection-free equilibrium
# ================================================================================================
DYN = ['T', 'EW', 'IW', 'EC', 'IC', 'P', 'VW', 'VT', 'C', 'H']       # the ten dynamic states (columns of rpt_adapt output)
def fig5():
    base = {k: v for k, v in read_input(os.path.join(HERE, 'input.txt')).items() if k not in ('w1', 'w2', 'w3', 'w4', 'C_crit', 'H_crit')}
    # (a) fixed-step RK4 (guard = 0), pure ODE (ext_floor = 0), same run as ./engine --selftest: TIP dose at tip_time, horizon t_end.
    #     Error = max over the ten states of |y_h - y_ref| / max(1, |y_ref|) at t_end, reference = a run with a much smaller step.
    hs_try = [0.05, 0.025, 0.0125, 0.00625, 0.003125]; h_ref = 0.0004
    def final(h):
        _, c = rpt_traj({**base, 'ext_floor': 0, 'guard': 0, 'h': h, 'dt_out': 0.25, 'adaptive': 0}); return [c[k][-1] for k in DYN], c['t'][-1]
    ref, tend = final(h_ref); errs = []
    for h in hs_try:
        y, _ = final(h); errs.append(max(abs(a - b) / max(1.0, abs(b)) for a, b in zip(y, ref)))
    orders = [math.log(errs[i] / errs[i + 1]) / math.log(2) for i in range(len(errs) - 1)]
    write_csv('fig05a_rk4_convergence', ['step_h_days', 'max_normalised_error'], list(zip(hs_try, errs)))
    # (b) return to the infection-free equilibrium: no virus (V_W0 = 0), no TIP, six perturbed starting states, 400 days.
    #     Deviation(t) = max over the ten states of |y(t) - y*| / max(1, |y*|); y* = the unperturbed run (checked to stay constant).
    eq_over = {**base, 'V_W0': 0, 'tip_dose': 0, 'ext_floor': 0, 'adaptive': 0, 't_end': 400, 'dt_out': 2}
    _, c0 = rpt_traj(eq_over); ystar = [c0[k][0] for k in DYN]
    drift = max(abs(c0[k][-1] - c0[k][0]) / max(1, abs(c0[k][0])) for k in DYN)
    pert = [('*T*(0) = 0.5 *T*^{*}', {'T0_factor': 0.5}, '#185FA5', None), ('*T*(0) = 2 *T*^{*}', {'T0_factor': 2}, '#185FA5', '8 5'),
            ('*C*(0) = 0', {'C0': 0}, '#D85A30', None), ('*C*(0) = 50 pg/mL', {'C0': 50}, '#D85A30', '8 5'),
            ('*H*(0) = 50 beats/min', {'H0_init': 50}, '#1D9E75', None), ('*H*(0) = 110 beats/min', {'H0_init': 110}, '#1D9E75', '8 5')]
    ser, rows = [], []
    for lab, o, col, dash in pert:
        _, c = rpt_traj({**eq_over, **o})
        dev = [max(abs(c[k][i] - ystar[j]) / max(1.0, abs(ystar[j])) for j, k in enumerate(DYN)) for i in range(len(c['t']))]
        ser.append(dict(x=c['t'], y=dev, color=col, dash=dash, label=lab)); rows += [(lab, t, d) for t, d in zip(c['t'], dev)]
    write_csv('fig05b_equilibrium_deviation', ['start_state', 'time_days', 'max_normalised_deviation'], rows)
    F = Fig('fig05', 1260, 640)
    F.text(14, 30, 'Numerical verification of the engine', S_TITLE, bold=True, kind='title')
    pa = [dict(x=hs_try, y=errs, color='#185FA5', width=2.6)]
    ref_line = dict(x=[hs_try[-1], hs_try[0]], y=[errs[-1] * 0.45, errs[-1] * 0.45 * (hs_try[0] / hs_try[-1]) ** 4], color='#6B6B6B', dash='8 5', width=1.8)
    px, py = panel(F, 8, 46, 620, 420, '(a) RK4 convergence against step size', dict(label='step size *h* (days)', log=True), dict(label='maximum normalised error', log=True),
                   [pa[0], ref_line], clamp=1e-300)
    for i, o in enumerate(orders):                  # local order printed under the data line at the segment mid-point
        xm = math.sqrt(hs_try[i] * hs_try[i + 1]); ym = math.sqrt(errs[i] * errs[i + 1])
        F.text(px(xm) - 12, py(ym) - 14, f'{o:.2f}', S_NOTE, 'end', fill='#185FA5', kind='note', note=True)
    px2, py2 = panel(F, 636, 46, 620, 420, '(b) Return to the infection-free equilibrium', dict(label='time (days)'), dict(label='max normalised deviation', log=True),
                     [dict(x=s['x'], y=s['y'], color=s['color'], dash=s['dash']) for s in ser], clamp=1e-12)
    F.text(636, 600, 'Deviations below 10^{−12} (numerical zero) are drawn on the 10^{−12} floor.', S_NOTE, fill='#333', kind='caption')
    y = legend(F, 14, 500, [('RK4, local order (labels beside the line)', '#185FA5', None), ('slope-4 reference', '#6B6B6B', '8 5')], 600)
    legend(F, 636, 500, [(s['label'], s['color'], s['dash']) for s in ser], 610)
    n = F.check(); F.save(os.path.join(OUT, 'fig05_numerical_verification.svg'))
    print(f'fig05: RK4 errors {["%.3e" % e for e in errs]} local orders {["%.3f" % o for o in orders]}; equilibrium drift without perturbation {drift:.2e}')
    return n

# ================================================================================================
# Figure 7: structure of the extended model (drawn from the right-hand sides in rpt_adapt.cpp; no numbers)
# ================================================================================================
def fig7():
    BW, BH, PITCH, X0 = 260, 72, 400, 110
    Lx = [X0 + PITCH * i for i in range(4)]; cx = [x + BW / 2 for x in Lx]
    FY0 = 52; yA0 = FY0 + 48; yB0 = yA0 + BH + 90; FY1 = yB0 + BH + 34; yC0 = FY1 + 130; yD0 = yC0 + BH + 100
    yA, yB, yC, yD = yA0 + BH / 2, yB0 + BH / 2, yC0 + BH / 2, yD0 + BH / 2
    F = Fig('fig07', 1740, yD0 + BH + 70)
    F.text(14, 32, 'Structure of the extended RPT-EBOV model', S_TITLE, bold=True, kind='title')
    F.frame(40, FY0, 1710, FY1); F.text(54, FY0 + 28, 'Infection and TIP dynamics (core)', S_LABEL, bold=True, kind='title')
    nm = {}
    def B(i, y, a, b, name): F.box(Lx[i], y, BW, BH, a, b, name); nm[name] = (Lx[i], y)
    B(0, yA0, '*T*', 'target cells', 'T');          B(1, yA0, '*E*_{W}', 'eclipse, WT-only', 'EW');          B(2, yA0, '*I*_{W}', 'productive, WT-only', 'IW');        B(3, yA0, '*V*_{W}', 'free WT virus', 'VW')
    B(0, yB0, '*P*', 'TIP-only primed', 'P');       B(1, yB0, '*E*_{Co}', 'eclipse, co-infected', 'ECo');    B(2, yB0, '*I*_{Co}', 'productive, co-infected', 'ICo');  B(3, yB0, '*V*_{T}', 'free TIP virus', 'VT')
    # row A / row B: forward flows (labels above the line, 12 px clear of it)
    def hlab(xa, xb, y, mk, above=True): F.text((xa + xb) / 2, y - 12 if above else y + 21, mk, S_LABEL, 'middle')
    F.arrow(40, yA, Lx[0], yA, name='lambda'); hlab(40, Lx[0], yA, '*λ*')
    for i, (rowy, labs) in enumerate([(yA, ['*β*_{W} *V*_{W}', '*k*_{W}']), (yB, ['*β*_{W} *V*_{W}', '*k*_{Co}'])]):
        for j in range(2):
            F.arrow(Lx[j] + BW, rowy, Lx[j + 1], rowy, name=f'row{i}{j}'); hlab(Lx[j] + BW, Lx[j + 1], rowy, labs[j])
    # TIP co-infection (vertical, beta_T V_T)
    for j, n in enumerate(['T', 'EW', 'IW']):
        F.arrow(cx[j], yA0 + BH, cx[j], yB0, name=f'tip{j}'); F.text(cx[j] + 12, (yA0 + BH + yB0) / 2 + 5, '*β*_{T} *V*_{T}', S_LABEL)
    # production: I_W -> V_W (p_W), I_Co -> V_T (p_T), and I_Co -> V_W (alpha p_W) joining the I_W -> V_W connector
    xj = Lx[2] + BW + 70
    F.arrow(Lx[2] + BW, yA, Lx[3], yA, name='pW'); F.text((Lx[2] + BW + Lx[3]) / 2 + 35, yA - 12, '*p*_{W}', S_LABEL, 'middle')
    F.arrow(Lx[2] + BW, yB, Lx[3], yB, name='pT'); F.text((Lx[2] + BW + Lx[3]) / 2, yB + 25, '*p*_{T}', S_LABEL, 'middle')
    F.arrow(xj, yB, xj, yA, head=None, name='alphapW', joins=('pW', 'pT')); F.text(xj + 12, (yA + yB) / 2 + 5, '*α p*_{W}', S_LABEL)
    # clearance (exits to the right)
    F.arrow(Lx[3] + BW, yA, Lx[3] + BW + 100, yA, name='cW'); F.text(Lx[3] + BW + 50, yA - 12, '*c*_{W}', S_LABEL, 'middle')
    F.arrow(Lx[3] + BW, yB, Lx[3] + BW + 100, yB, name='cT'); F.text(Lx[3] + BW + 50, yB - 12, '*c*_{T}', S_LABEL, 'middle')
    # bottom row: cytokine (col 0), adaptive block B, A, Z (cols 1..3); pulse under the cytokine
    B(0, yC0, '*C*', 'cytokine', 'C'); B(0, yD0, '*H*', 'pulse rate', 'H')
    B(1, yC0, '*B*', 'antibody-secreting pool', 'B'); B(2, yC0, '*A*', 'neutralising antibody', 'A'); B(3, yC0, '*Z*', 'cytotoxic T-cell pool', 'Z')
    ymid = (FY1 + yC0) / 2 + 5
    # core <-> cytokine: production down (solid), cytokine-enhanced death of infected cells up (dashed = indirect effect)
    F.arrow(cx[0] - 50, FY1, cx[0] - 50, yC0, name='rho'); F.text(cx[0] - 62, ymid, '*ρ*_{I}, *ρ*_{V}, *ρ*_{T}', S_LABEL, 'end')
    F.arrow(cx[0] + 50, yC0, cx[0] + 50, FY1, dash=True, name='dN'); F.text(cx[0] + 62, ymid, '*δ*_{N}', S_LABEL)
    # core -> adaptive (drive), adaptive -> core (flat end = inhibition)
    F.arrow(cx[1], FY1, cx[1], yC0, name='piB'); F.text(cx[1] + 12, ymid, '*π*_{B} (*V*_{W})', S_LABEL)
    F.arrow(cx[2] - 50, yC0, cx[2] - 50, FY1, head='flat', name='epsAb'); F.text(cx[2] - 62, ymid, '*ε*_{Ab}', S_LABEL, 'end')
    F.arrow(cx[2] + 50, yC0, cx[2] + 50, FY1, head='flat', name='epsClr'); F.text(cx[2] + 62, ymid, '*ε*_{Clr}', S_LABEL)
    F.arrow(cx[3] - 50, FY1, cx[3] - 50, yC0, name='piZ'); F.text(cx[3] - 62, ymid, '*π*_{Z} (*I*)', S_LABEL, 'end')
    F.arrow(cx[3] + 50, yC0, cx[3] + 50, FY1, head='flat', name='muZ'); F.text(cx[3] + 62, ymid, '*μ*_{Z}', S_LABEL)
    F.arrow(Lx[1] + BW, yC, Lx[2], yC, name='pA'); F.text((Lx[1] + BW + Lx[2]) / 2, yC - 12, '*p*_{A}', S_LABEL, 'middle')
    # cytokine <-> pulse (vagal feedback): drive down (solid), vagal gain m(H) back up (dashed)
    ym2 = (yC0 + BH + yD0) / 2 + 5
    F.arrow(cx[0] - 30, yC0 + BH, cx[0] - 30, yD0, name='dH'); F.text(cx[0] - 42, ym2, 'Δ_{H}', S_LABEL, 'end')
    F.arrow(cx[0] + 30, yD0, cx[0] + 30, yC0 + BH, dash=True, name='mH'); F.text(cx[0] + 42, ym2, '*m*(*H*)', S_LABEL)
    # group captions and key
    F.text(cx[2], yC0 + BH + 32, 'Adaptive block (Section 11)', S_LABEL, 'middle', bold=True, kind='caption')
    F.text(Lx[0], yD0 + BH + 32, 'Cytokine and pulse', S_LABEL, 'start', bold=True, kind='caption')
    kx, ky = Lx[1], yD0 + 24
    F.arrow(kx, ky, kx + 50, ky, name='key1'); F.text(kx + 64, ky + 5, 'conversion, production or drive', S_LABEL, kind='key')
    F.arrow(kx + 400, ky, kx + 450, ky, head='flat', name='key2'); F.text(kx + 464, ky + 5, 'inhibition (flat end)', S_LABEL, kind='key')
    F.arrow(kx + 700, ky, kx + 750, ky, dash=True, name='key3'); F.text(kx + 764, ky + 5, 'indirect effect (dashed)', S_LABEL, kind='key')
    n = F.check(); F.save(os.path.join(OUT, 'fig07_model_schematic.svg')); return n

# ================================================================================================
# Figure 13: naive and vaccinated hosts, untreated and with the TIP profile (A_vac = 50 ug/mL, c_T = 0.5 /d)
# ================================================================================================
def fig13():
    prof = {k: v for k, v in AS.BEST.items()}; AVAC = 50.0; CT = 0.5
    runs = [('naive, untreated', {'tip_dose': 0}, '#6B6B6B', None), ('naive + TIP', {'c_T': CT}, '#D85A30', '8 5'),
            ('vaccinated, untreated', {'tip_dose': 0, 'A0': AVAC}, '#185FA5', None), ('vaccinated + TIP', {'c_T': CT, 'A0': AVAC}, '#1D9E75', '8 5')]
    res = []
    for lab, o, col, dash in runs:
        s, c = rpt_traj({**prof, **o, 'adaptive': 1}); res.append((lab, col, dash, s, c))
    rows = [(lab, t, v, a, cy, b) for lab, _, _, _, c in res for t, v, a, cy, b in zip(c['t'], c['VW'], c['A'], c['C'], c['B'])]
    write_csv('fig13_trajectories', ['run', 'time_days', 'V_W', 'A', 'C', 'B'], rows)
    F = Fig('fig13', 1260, 880); F.text(14, 30, 'Naive and vaccinated hosts, untreated and with the TIP profile', S_TITLE, bold=True, kind='title')
    cfg = [('(a) Free WT virus (log scale)', 'VW', 'copies/mL', True, [(2e-4, '*V*_{ext}', '#A32D2D', '2 4')]), ('(b) Neutralising antibody', 'A', 'antibody *A* (µg/mL)', False, []),
           ('(c) Cytokine', 'C', 'cytokine (pg/mL)', False, [(500, '*C*_{crit}', '#A32D2D', '2 4')]), ('(d) Antibody-secreting pool', 'B', '*B* (relative units)', False, [])]
    for k, (ttl, key, yl, lg, hl) in enumerate(cfg):
        X, Y = 8 + (k % 2) * 628, 46 + (k // 2) * 360
        panel(F, X, Y, 620, 350, ttl, dict(label='time (days)'), dict(label=yl, log=lg), [dict(x=c['t'], y=c[key], color=col, dash=dash) for _, col, dash, _, c in res], hlines=hl)
    legend(F, 14, 800, [(lab, col, dash) for lab, col, dash, _, _ in res], 1200)
    n = F.check(); F.save(os.path.join(OUT, 'fig13_vaccination_trajectories.svg'))
    for lab, _, _, s, c in res: print(f"fig13 {lab:24s} peakV={s['peakV']:.3g} peakC={s['peakC']:.3g} extinct={int(s['extinct'])} tExt={s['tExt']:.3g} A(30)={c['A'][-1]:.3g}")
    return n

# ================================================================================================
# Figure 15: peak virus and peak cytokine against vaccine-induced antibody A_vac
# ================================================================================================
def fig15():
    prof = dict(AS.BEST); grid = [2.5 * i for i in range(0, 41)]     # A_vac = 0 ... 100 ug/mL in steps of 2.5
    arms = [('untreated', {'tip_dose': 0}, '#6B6B6B', None), ('TIP, *c*_{T} = 0.5 d^{−1}', {'c_T': 0.5}, '#D85A30', None), ('TIP, *c*_{T} = 0.1 d^{−1}', {'c_T': 0.1}, '#1D9E75', '8 5')]
    ser, rows = [], []
    for lab, o, col, dash in arms:
        pv, pc = [], []
        for a in grid:
            s = AS.run({**prof, **o, 'A0': a, 'adaptive': 1}); pv.append(s['peakV']); pc.append(s['peakC']); rows.append((lab, a, s['peakV'], s['peakC']))
        ser.append((lab, col, dash, pv, pc))
    write_csv('fig15_antibody_sweep', ['arm', 'A_vac_ug_per_mL', 'peak_V_W', 'peak_C'], rows)
    F = Fig('fig15', 1260, 520); F.text(14, 30, 'Vaccine-induced antibody: peak virus and peak cytokine', S_TITLE, bold=True, kind='title')
    panel(F, 8, 46, 620, 400, '(a) Peak free WT virus', dict(label='vaccine-induced antibody *A*_{vac} (µg/mL)', lo=0, hi=100), dict(label='peak *V*_{W} (copies/mL)', log=True),
          [dict(x=grid, y=pv, color=col, dash=dash) for _, col, dash, pv, _ in ser], hlines=[(AS.BEST and 100.0, 'inoculum', '#A32D2D', '2 4')])
    panel(F, 636, 46, 620, 400, '(b) Peak cytokine', dict(label='vaccine-induced antibody *A*_{vac} (µg/mL)', lo=0, hi=100), dict(label='peak cytokine (pg/mL)'),
          [dict(x=grid, y=pc, color=col, dash=dash) for _, col, dash, _, pc in ser], hlines=[(500, '*C*_{crit}', '#A32D2D', '2 4')])
    legend(F, 14, 482, [(lab, col, dash) for lab, col, dash, _, _ in ser], 1200)
    n = F.check(); F.save(os.path.join(OUT, 'fig15_vaccine_antibody_sweep.svg'))
    for lab, _, _, pv, pc in ser: print(f'fig15 {lab:28s} peakV(A=0)={pv[0]:.3g} peakV(A=60)={pv[grid.index(60.0)]:.3g} peakC(A=0)={pc[0]:.3g} peakC(A=60)={pc[grid.index(60.0)]:.3g}')
    return n

if __name__ == '__main__':
    bad = sum(f() for f in (fig5, fig7, fig13, fig15))
    print('make_figures: all layout checks passed' if bad == 0 else f'make_figures: {bad} layout problem(s)')
    sys.exit(1 if bad else 0)
