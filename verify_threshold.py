#!/usr/bin/env python3
"""Threshold check behind the 'seven digits' statement (Sections 7.5 and 15).
Bisects on the sign of the largest eigenvalue of the (E_W, I_W, V_W) block (Eq 16) for the entry-blocking
effect and for the clearance gain, and compares the crossing points with the analytic values 1 - 1/R0
and R0 - 1 (Eq 12-13). Reads input.txt; needs numpy only. Prints 9 digits so the agreement can be read off.
Usage: python3 verify_threshold.py [input.txt]"""
import sys, re
import numpy as np
f = sys.argv[1] if len(sys.argv) > 1 else 'input.txt'
P = {}
for line in open(f):
    line = line.split('#')[0]
    if '=' in line:
        k, v = line.split('=', 1)
        try: P[k.strip()] = float(v)
        except ValueError: pass
Ts = P['lambda'] / P['d_T']; Cs = P['C_rest']
cyt = P['delta_N'] * Cs / (P['K_F'] + Cs)
d = P['delta_W'] + cyt
R0 = P['beta_W'] * P['p_W'] * Ts / (P['c_W'] * d)                     # Eq 12
def lam(eAb, eCl):                                                     # largest real part of the linearised block
    bt = (1 - eAb) * P['beta_W'] * Ts; k = P['k_W']; c = P['c_W'] * (1 + P['kappa'] * eCl)
    return max(np.linalg.eigvals(np.array([[-k, 0, bt], [k, -d, 0], [0, P['p_W'], -c]])).real)
def bisect(g, lo, hi):
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if g(mid) > 0: lo = mid                                        # still growing: threshold is higher
        else: hi = mid
    return 0.5 * (lo + hi)
a = bisect(lambda e: lam(e, 0), 0.0, 0.999999)
b = bisect(lambda e: lam(0, e), 0.0, 1e3)
print(f'R0 = {R0:.9f}')
print(f'entry blocking : bisection {a:.9f}   analytic 1-1/R0 = {1 - 1 / R0:.9f}   |diff| = {abs(a - (1 - 1 / R0)):.2e}')
print(f'clearance gain : bisection {b:.9f}   analytic R0-1   = {R0 - 1:.9f}   |diff| = {abs(b - (R0 - 1)):.2e}')
ok = abs(a - (1 - 1 / R0)) < 5e-8 and abs(b - (R0 - 1)) < 5e-7
print('threshold match to seven significant digits:', 'PASS' if ok else 'FAIL')
sys.exit(0 if ok else 1)
