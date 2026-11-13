#!/usr/bin/env python3
"""Put the regenerated figures and reproducible captions into the manuscript .docx.
Usage: python3 build_docx_figures.py [docx] [fig numbers ...]   (default: all 17)
Needs rsvg-convert. Only caption text and picture files are touched."""
import re, sys, zipfile, subprocess, os, shutil, tempfile
HERE = os.path.dirname(os.path.abspath(__file__))
DOCX = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, 'Rein_Persistence_Theory_EBOV_Nyachio_Final.docx')
ONLY = {int(a) for a in sys.argv[2:]} or set(range(1, 18))
SVG = {1:'sample_best.svg',2:'fig02_extinction_time_vs_cT.svg',3:'screen_pareto.svg',4:'fig04_weak_vs_best.svg',
 5:'fig05_numerical_verification.svg',6:'fig06_threshold_growth_rate.svg',7:'fig07_model_schematic.svg',
 8:'fig08_extinction_map.svg',9:'fig09_trajectories_near_threshold.svg',10:'fig10_threshold_vs_design.svg',
 11:'fig11_local_sensitivity.svg',12:'fig12_cytokine_tradeoff.svg',13:'fig13_vaccination_trajectories.svg',
 14:'fig14_vaccination_threshold.svg',15:'fig15_vaccine_antibody_sweep.svg',16:'fig16_establishment_probability.svg',
 17:'fig17_stochastic_tail.svg'}
# Captions: text after "Figure N." ; markup _{sub} ^{sup}. Every number is printed by the scripts (verification_log.txt, figure_data/).
R = lambda f, s='make_figures_analysis.py': f' Reproduce: {s} {f}; data in figure_data/.'
CAP = {
2:' Time to WT extinction against free-TIP clearance rate c_{T} for three TIP profiles. Strongest: refined profile (analysis_scripts.REFINED) at dose 10^{9}, day 1; Best: the screened best profile (screen_report.txt); Weaker: the WEAK profile of make_figures_analysis.py (β_{T} = 5×10^{−7}, α = 0.2, p_{T} = 1500, dose 10^{8}, day 1). Solid: primed cells tracked (reference); dashed: abortive-entry variant; triangles: thresholds c_{T}^{*} (primed 0.456 and 0.256; abortive 1.151, 0.696 and 0.525; the weaker profile has no primed extinction); shaded: WT clearance range 5–15 d^{−1}.'+R('fig2'),
4:' Untreated infection against a weaker TIP profile (WEAK in make_figures_analysis.py: β_{T} = 5×10^{−7}, α = 0.2, p_{T} = 1500, dose 10^{8}, day 1) and the best screened profile (c_{T} = 0.18 d^{−1}; screen_report.txt), illustrative parameters: (a) WT virus over time; (b) peak load (untreated 1.3×10^{7}, weaker 4.7×10^{5}, best 100 = inoculum, i.e. extinct); (c) peak cytokine against its safety bound (207, 31 and 5.0 against C_{crit} = 500).'+R('fig4'),
5:' Numerical verification. (a) Maximum error of RK4 against step size h = 0.05 to 0.003125 d (reference h = 0.0004 d, 30 days, extinction floor and guard off), with the slope-4 reference; local orders 4.33, 4.17, 4.09 and 4.05 between successive steps. (b) Maximum normalised deviation from the infection-free equilibrium over 400 days for six perturbed starting states (T_{0} factor 0.5 or 2, C_{0} 0 or 50, H_{0} 50 or 110), no virus and no TIP; deviations below 10^{−12} are drawn on the floor.'+R('fig5','make_figures.py'),
8:' Extinction map in (c_{T}, dose) space for the fixed profile, a 36 × 48 grid (c_{T} 0.01–20 d^{−1}, dose 10^{5}–10^{9}). Green: WT extinct within 30 days; beige: no extinction. (a) Primed cells tracked (reference model); (b) abortive-entry variant. The grey band is the WT clearance range, 5–15 d^{−1}.'+R('fig8'),
9:' Trajectories for seven values of c_{T} (0.05, 0.1, 0.3, 0.5, 1, 2 and 5 d^{−1}; colour key below the panels) with the fixed profile and the extinction floor on. (a) Free WT virus with the extinction threshold V_{ext}; (b) free TIP; (c) phase portrait of V_{W} against V_{T}. Extinction occurs for c_{T} = 0.05 and 0.1 (day 13.3); all larger values regrow, consistent with c_{T}^{*} = 0.256 d^{−1}.'+R('fig9'),
10:' Persistence threshold c_{T}^{*} against (a) hijacking factor α, (b) production ratio p_{T}/p_{W}, (c) entry ratio β_{T}/β_{W} and (d) TIP dose, with primed cells tracked (solid) and in the abortive variant (dashed). Grey band: WT clearance range. In (a) the reference model has no extinction for α of 0.095 and above (last extinct grid point α = 0.07).'+R('fig10'),
12:' Trade-off between persistence threshold (blue, left axis) and peak cytokine at the threshold (red, right axis) as TIP innate stimulation ρ_{T} increases over 13 values from 0 to 3×10^{−3}; c_{T}^{*} rises from 0.256 to 0.373 d^{−1} and peak C from 5 to 1.7×10^{4}. The dotted line is the cytokine bound C_{crit}. The leftmost point is ρ_{T} = 0.'+R('fig12'),
14:' Vaccination and the persistence threshold. (a) c_{T}^{*} against vaccine-induced antibody for four assumptions on antibody recognition of the TIP (ξ = 0, 0.5, 1 multiplicative and 1 additive); the grey band is the WT clearance range and the dotted line the naive-host threshold. (b) R_{eff}(A) for vaccine alone (Eq 13); dotted lines: R_{eff} = 1 at A = 68.6 µg mL^{−1} and 60-day extinction at 94.0 µg mL^{−1}.'+R('fig14'),
16:' Establishment of infection. (a) Monte Carlo fraction of lineages that establish (±2 standard errors) against the prediction 1 − q of Eq 26 for the eight parameter sets of Table 16; the dashed line is equality. (b) Probability of establishing infection against inoculum size for a naive host (R_{eff} = 8.73; P = 0.27 for 10 virions), a vaccinated host (A = 50 µg mL^{−1}, R_{eff} = 1.19; P = 0.008), R_{eff} = 3 and R_{eff} = 0.9.'+R('fig16'),
17:' Probability of WT extinction in the stochastic tail against free-TIP clearance, at 11 values of c_{T}, with 95% Wilson confidence intervals (200 replicates per point, 300 starting virions), and the deterministic outcome (dashed). The deterministic threshold c_{T}^{*} = 0.256 d^{−1} is marked.'+R('fig17'),
}
def esc(s): return s.replace('&','&amp;').replace('<','&lt;').replace('>','&gt;')
def runs(text, rpr_base, extra=''):
    out=[]
    for m in re.finditer(r'([_^])\{([^}]*)\}|([^_^]+|[_^])', text):
        if m.group(1):
            va='subscript' if m.group(1)=='_' else 'superscript'
            out.append(f'<w:r><w:rPr>{rpr_base}<w:vertAlign w:val="{va}"/></w:rPr><w:t xml:space="preserve">{esc(m.group(2))}</w:t></w:r>')
        else:
            out.append(f'<w:r><w:rPr>{rpr_base}</w:rPr><w:t xml:space="preserve">{esc(m.group(3))}</w:t></w:r>')
    return ''.join(out)
tmp = tempfile.mkdtemp()
zin = zipfile.ZipFile(DOCX)
doc = zin.read('word/document.xml').decode('utf8')
rels = zin.read('word/_rels/document.xml.rels').decode('utf8')
pm = list(re.finditer(r'<w:p[ >].*?</w:p>', doc, flags=re.S))
def ptext(p): return ''.join(re.findall(r'<w:t[^>]*>([^<]*)</w:t>', p))
cap = {}
for i, m in enumerate(pm):
    mm = re.match(r'Figure (\d+)\.', ptext(m.group(0)))
    if mm: cap[int(mm.group(1))] = i
newmedia = {}; repl = []
for n in sorted(ONLY):
    ci = cap[n]; cp = pm[ci].group(0); dp = pm[ci-1].group(0)
    if n in CAP:
        runs_ = re.findall(r'<w:r>.*?</w:r>', cp, flags=re.S)
        pre = cp[:cp.index('<w:r>')]
        bold = '<w:b/>' in runs_[0] and len(runs_) > 1
        if bold:
            rb = re.search(r'<w:rPr>(.*?)</w:rPr>', runs_[1], flags=re.S).group(1)
            rb = re.sub(r'<w:vertAlign[^>]*/>', '', rb)
            body = runs_[0] + runs(CAP[n], rb)
        else:
            rb = re.sub(r'<w:vertAlign[^>]*/>', '', re.search(r'<w:rPr>(.*?)</w:rPr>', runs_[0], flags=re.S).group(1))
            body = runs(f'Figure {n}.' + CAP[n], rb)
        repl.append((cp, pre + body + '</w:p>'))
    rid = re.search(r'r:embed="(rId\d+)"', dp)
    if not rid: print('no picture before caption', n); continue
    tgt = re.search(rf'Id="{rid.group(1)}"[^>]*Target="([^"]+)"|Target="([^"]+)"[^>]*Id="{rid.group(1)}"', rels)
    tgt = tgt.group(1) or tgt.group(2)
    png = os.path.join(tmp, f'f{n}.png')
    subprocess.check_call(['rsvg-convert', '-z', '2', os.path.join(HERE, SVG[n]), '-o', png])
    from struct import unpack
    b = open(png, 'rb').read(); w, h = unpack('>II', b[16:24])
    newmedia['word/' + tgt] = b
    cx = int(re.search(r'<wp:extent cx="(\d+)"', dp).group(1)); cy = round(cx * h / w)
    nd = re.sub(r'(<wp:extent cx="\d+" cy=")\d+', rf'\g<1>{cy}', dp)
    nd = re.sub(r'(<a:ext cx="\d+" cy=")\d+', rf'\g<1>{cy}', nd)
    repl.append((dp, nd))
for a, b in repl: doc = doc.replace(a, b, 1)
out = DOCX + '.new'
with zipfile.ZipFile(out, 'w', zipfile.ZIP_DEFLATED) as zo:
    for it in zin.infolist():
        data = zin.read(it.filename)
        if it.filename == 'word/document.xml': data = doc.encode('utf8')
        elif it.filename in newmedia: data = newmedia[it.filename]
        zo.writestr(it, data)
os.replace(out, DOCX); shutil.rmtree(tmp)
print('figures written', sorted(ONLY))
