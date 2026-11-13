// =====================================================================================
// code.cpp  -  EBOV-TIP within-host engine and TIP-profile screening  (C++20, std library only)
// Reference implementation of the Rein Persistence Theory of Persistence-Limited
// Interference in Orthoebolavirus zairense (Ebola virus, EBOV).      Author: Reinford Nyachio
//
// WHAT THIS FILE IS
//   A deterministic ODE simulator of one infected host: wild-type (WT) Ebola virus plus an
//   optional dose of Therapeutic Interfering Particles (TIPs), a cytokine surrogate and a pulse
//   rate. On top of the simulator sits a screening module that tries thousands of candidate
//   TIP "profiles" and ranks them by efficacy and safety. Equations 1-14 of the proposal are
//   implemented literally; each line of rhs() carries its equation number.
//
// HOW TO BUILD AND RUN (any C++20 compiler; no third-party libraries)
//   g++ -std=c++20 -O2 -pthread -Wall -Wextra code.cpp -o engine
//   ./engine [input.txt] [output.txt] [errors.txt]            one trajectory (+ output.svg)
//   ./engine --screen [screen.txt] [input.txt] [errors.txt]   TIP-profile screening
//   ./engine --selftest [input.txt]                           RK4 convergence-order check
//   g++ -std=c++20 -O2 -pthread -DEBOV_NO_MAIN test_api.cpp code.cpp -o test_api   (C-API test)
//
// STATE VECTOR: 10 dynamic states (T, EW, IW, VW, EC, IC, VT, C, H, P) followed by 4 running
//   integrals (JV, JC, JH, U) that accumulate the terms of the objective J (Eq 14), so the
//   objective weights can be re-applied after a run without re-simulating.
//
// UNITS: cells/mL, virions (copies)/mL, time in DAYS, cytokine in pg/mL, pulse in beats/min.
//
// READING ORDER FOR A NEWCOMER
//   1 Params, rhs, rk4        the mathematics      2 run()        one simulation
//   3 evaluate, screen()      the search           4 C API        how Kotlin/JNI calls it
//
// STATUS NOTES (what the proposal plans but this file does NOT yet contain)
//   * No 64-byte alignment / structure-of-arrays tiling: State is a plain std::array<double,14>.
//   * No JNI glue: only the C API at the bottom; the Kotlin/JNI layer is a Phase-3 deliverable.
//   * R_eff here is the entry+clearance value only; the next-generation-matrix R_eff with a TIP
//     dose (proposal Section 4.4) and the Gillespie stochastic extinction check are not coded.
//   * Workers use std::thread, not std::jthread (works on every NDK toolchain).
//   * Every default parameter is an ILLUSTRATIVE placeholder, not a calibrated measurement.
//   * Parameter structs are read/written by byte offset in evaluate()/load_screen_stream():
//     Params must stay a plain struct whose swept members are doubles (see notes there).
// =====================================================================================
// Only the standard library is used (zero-dependency commitment): a future reader needs a
// compiler and nothing else. <thread>/<atomic> provide the parallel screening workers.
#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <limits>
#include <map>
#include <numeric>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

// Index names for the state vector. Keeping them as an enum means every equation below reads
// like the paper: y[iVW] is V_W, y[iC] is the cytokine C, and so on. N is the vector length.
//   T  uninfected target cells        EW/IW  eclipse/productive WT-only infected cells
//   EC/IC  eclipse/productive co-infected cells (WT + TIP)      P  TIP-only "primed" cells
//   VW/VT  free WT / free TIP virions    C  cytokine (IL-6-like)    H  pulse rate
// 0..9 dynamic states (9 = P, TIP-only primed cells), 10..13 running integrals (JV, JC, JH, U) so objective weights can be applied afterwards
enum Idx { iT, iEW, iIW, iVW, iEC, iIC, iVT, iC, iH, iP, iJV, iJC, iJH, iU, N };
// kDyn = number of true dynamic states (the first ten); the last four entries are bookkeeping integrals.
constexpr int kDyn = 10;
// One point in time. A plain array: copying it is cheap and contains no hidden state.
using State = std::array<double, N>;

// Every number the model needs, with its default (the illustrative reference set; input.txt
// overrides any of them). Naming: b = beta (infection rate), d = delta/death, k = eclipse
// exit rate, p = production, c = clearance. W = wild type, T = TIP, Co = co-infected.
// Line 1: cell kinetics.  Line 2: death/production/hijacking/clearance.  Line 3: cytokine-
// mediated killing and the antibody-like entry (eps_ab) and clearance (eps_clr) effects.
// Line 4-5: cytokine and autonomic (vagal/pulse) subsystem.  Line 6: objective weights and
// safety bounds.  Line 7: run control (initial inoculum, TIP dose/time, horizon, step).
// The last line holds values DERIVED from the others (see derive()); never set them by hand.
struct Params {
  double lambda = 5e4, dT = 0.07, bW = 5e-7, bT = 5e-7, kW = 6, kCo = 6;
  double dW = 2, dCo = 2, pW = 500, pT = 1500, alpha = 0.2, cW = 10, cT = 10;
  double dN = 0.5, KF = 50, eps_ab = 0, eps_clr = 0, kappa = 1;
  double c_rest = 5, dC = 5, gamma = 3, KC = 50, rhoI = 0.005, rhoV = 1e-6;
  double sigma_v = 0.5, H0 = 70, Hmax = 130, dH = 50, KH = 100, n_hill = 2, tau_H = 0.02;
  double w1 = 1, w2 = 1, w3 = 1, w4 = 1e-6, C_crit = 500, H_crit = 130, v_ext = 2e-4;
  double primed = 1, rhoT = 0, ext_floor = 0, VW0 = 100, tip_dose = 0, tip_time = 0, t_end = 30, h = 0.005, dt_out = 0.25;
  double Ts = 0, Cs = 0, Hs = 0, s0 = 0;  // derived
};

// Hill function x^n/(K^n+x^n): a smooth saturating switch from 0 to 1 centred on K.
// Integer n is used so the power is exact multiplication (fast, deterministic, no pow()).
static double hill(double x, double K, int n) {
  double a = 1, b = 1;
  for (int i = 0; i < n; ++i) { a *= x; b *= K; }
  return a / (b + a);
}
// Vagal (cholinergic anti-inflammatory) gain m(H) of Eq 9: equal to 1-sigma_v at resting pulse
// H0 (maximum suppression of cytokine synthesis) and fading to 1 as H rises to Hmax.
// Bounded in [1-sigma_v, 1] by construction, which keeps the feedback loop well behaved.
static double vagal(const Params& p, double H) {
  return 1.0 - p.sigma_v * std::clamp((p.Hmax - H) / (p.Hmax - p.H0), 0.0, 1.0);
}

// THE MODEL. Right-hand side dy/dt = f(y) of Equations 1-10 (+14). PURE: no globals, no I/O,
// no randomness, so the same input always gives the same output bit for bit.
//   cyt  = cytokine-mediated extra death rate of infected cells (saturating in C)
//   inf  = rate of new WT infection of uninfected cells          (entry blocked by 1-eps_ab)
//   tip  = per-cell rate at which free TIP enters a cell          (beta_T * V_T)
//   supW = per-cell rate at which free WT enters a cell (used for superinfection of P cells)
static void rhs(const Params& p, const State& y, State& d) {
  // ADVANCED EXTENSIONS (see ADVANCED_EXTENSIONS.md & PDF):
  // Multi-stage entry would insert intermediate states B (bound) and K (NPC1-engaged)
  // before E_W; β_eff would multiply (1-ε_drug). Primed superinfection gains factor (1-π).
  // Cytokine rhs gains / (1+R/K_R) and optional g_cond term. Latent L seeds from E_W+I_W.
  // None of these states are present in the current 10-state vector; setting the advanced
  // parameters to their null values recovers this rhs exactly.
  const double ent = 1.0 - p.eps_ab;
  const double cyt = p.dN * y[iC] / (p.KF + y[iC]);
  const double inf = ent * p.bW * y[iT] * y[iVW];
  const double tip = ent * p.bT * y[iVT];
  const double supW = ent * p.bW * y[iVW];                                       // WT entry into a primed cell
  // --- cells (Eq 1-5, 4b). TIP entering an uninfected cell makes a PRIMED cell P (it cannot
  // replicate the TIP alone); WT can later enter P and turn it into a co-infected eclipse cell.
  // That superinfection path is the heart of the primed-cell postulate; primed=0 removes it.
  d[iT]  = p.lambda - p.dT * y[iT] - inf - tip * y[iT];                          // Eq 1
  d[iEW] = inf - p.kW * y[iEW] - tip * y[iEW];                                   // Eq 2
  d[iIW] = p.kW * y[iEW] - (p.dW + cyt) * y[iIW] - tip * y[iIW];                 // Eq 3
  d[iEC] = tip * y[iEW] + p.primed * supW * y[iP] - p.kCo * y[iEC];              // Eq 4 (+ superinfection of primed cells)
  d[iIC] = p.kCo * y[iEC] + tip * y[iIW] - (p.dCo + cyt) * y[iIC];               // Eq 5
  d[iP]  = p.primed * (tip * y[iT] - supW * y[iP]) - p.dT * y[iP];               // Eq 4b primed (TIP-only) cells
  // --- virus (Eq 6-7). Inside co-infected cells WT output is cut to alpha*pW (hijacking) while
  // TIP output is pT (> pW). Clearance of WT is boosted by antibody-like eps_clr; free TIP is
  // cleared at cT. The small ratio cT/cW is what the Persistence Theory is about.
  d[iVW] = p.pW * y[iIW] + p.alpha * p.pW * y[iIC]
         - p.cW * (1.0 + p.kappa * p.eps_clr) * y[iVW];                          // Eq 6
  d[iVT] = p.pT * y[iIC] - p.cT * y[iVT];                                        // Eq 7
  // --- cytokine (Eq 8-9) and pulse (Eq 10). m(H) multiplies ALL synthesis terms: basal s0,
  // infected-cell term rhoI, virus term rhoV, optional TIP-innate term rhoT (default 0) and the
  // second-order Hill self-amplification gamma*C^2/(KC^2+C^2); first-order clearance dC.
  d[iC]  = vagal(p, y[iH]) * (p.s0 + p.rhoI * (y[iIW] + y[iIC]) + p.rhoV * y[iVW] + p.rhoT * y[iVT]
         + p.gamma * hill(y[iC], p.KC, 2)) - p.dC * y[iC];                       // Eq 8-9
  d[iH]  = (p.H0 + p.dH * hill(y[iC], p.KH, (int)p.n_hill) - y[iH]) / p.tau_H;  // Eq 10
  // --- running integrals for the objective J (Eq 14): deviations are scaled by the safety
  // bounds and taken relative to the infection-free equilibrium (Cs, Hs); U = TIP exposure.
  const double ec = (y[iC] - p.Cs) / p.C_crit, eh = (y[iH] - p.Hs) / p.H_crit;
  d[iJV] = std::log10(1.0 + std::max(y[iVW], 0.0));                              // Eq 14 terms
  d[iJC] = ec * ec;
  d[iJH] = eh * eh;
  d[iU]  = y[iVT];
}

// Classical 4th-order Runge-Kutta, one fixed step h. Fixed (not adaptive) on purpose: with the
// same step sequence, results are reproducible bit for bit on any machine/thread count. The
// step is validated by the step-halving self-test (order should be ~4) and by jac_bound().
static void rk4(const Params& p, State& y, double h) {
  State k1, k2, k3, k4, s;
  rhs(p, y, k1);
  for (int i = 0; i < N; ++i) s[i] = y[i] + 0.5 * h * k1[i];
  rhs(p, s, k2);
  for (int i = 0; i < N; ++i) s[i] = y[i] + 0.5 * h * k2[i];
  rhs(p, s, k3);
  for (int i = 0; i < N; ++i) s[i] = y[i] + h * k3[i];
  rhs(p, s, k4);
  for (int i = 0; i < N; ++i) y[i] += (h / 6.0) * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]);
}

// STABILITY GUARD. RK4 is stable only if h*|lambda_max| < ~2.78 (real axis). Infection terms
// such as bW*VW can reach tens-hundreds per day at peak viraemia, so before each integration
// segment we bound the Jacobian spectral radius by its infinity norm (rows of |df_i/dy_j|,
// by forward differences) and advance() shrinks h if h*bound would exceed 2.5.
// infinity-norm of a finite-difference Jacobian of the dynamic block: bound on the spectral radius
static double jac_bound(const Params& p, const State& y) {
  State f0, f1, z = y;
  rhs(p, y, f0);
  std::array<double, kDyn> row{};
  for (int j = 0; j < kDyn; ++j) {
    const double e = 1e-6 * std::max(1.0, std::fabs(y[j]));
    z[j] = y[j] + e;
    rhs(p, z, f1);
    z[j] = y[j];
    for (int i = 0; i < kDyn; ++i) row[i] += std::fabs((f1[i] - f0[i]) / e);
  }
  return *std::max_element(row.begin(), row.end());
}

// ---- Parameter input. Plain-text "name = value" lines; '#' starts a comment. Unknown names,
// bad numbers and out-of-range values are REPORTED, never silently ignored.
// ------------------------------------------------------------------ parameters
// Table row linking a text-file name (e.g. "beta_T") to the C++ member it sets.
struct Entry { const char* name; double* ptr; };
// The single source of truth for which parameters can be set from text and swept in screening.
// To add a parameter: add the member to Params, one row here, and (if needed) a check in validate().
static std::vector<Entry> registry(Params& p) {
  return {{"lambda", &p.lambda}, {"d_T", &p.dT}, {"beta_W", &p.bW}, {"beta_T", &p.bT},
    {"k_W", &p.kW}, {"k_Co", &p.kCo}, {"delta_W", &p.dW}, {"delta_Co", &p.dCo},
    {"p_W", &p.pW}, {"p_T", &p.pT}, {"alpha", &p.alpha}, {"c_W", &p.cW}, {"c_T", &p.cT},
    {"delta_N", &p.dN}, {"K_F", &p.KF}, {"eps_Ab", &p.eps_ab}, {"eps_Clr", &p.eps_clr},
    {"kappa", &p.kappa}, {"C_rest", &p.c_rest}, {"d_C", &p.dC}, {"gamma", &p.gamma},
    {"K_C", &p.KC}, {"rho_I", &p.rhoI}, {"rho_V", &p.rhoV}, {"rho_T", &p.rhoT}, {"primed", &p.primed}, {"sigma_v", &p.sigma_v},
    {"H_0", &p.H0}, {"H_max", &p.Hmax}, {"Delta_H", &p.dH}, {"K_H", &p.KH},
    {"n_hill", &p.n_hill}, {"tau_H", &p.tau_H}, {"w1", &p.w1}, {"w2", &p.w2},
    {"w3", &p.w3}, {"w4", &p.w4}, {"C_crit", &p.C_crit}, {"H_crit", &p.H_crit},
    {"V_ext", &p.v_ext}, {"ext_floor", &p.ext_floor}, {"V_W0", &p.VW0}, {"tip_dose", &p.tip_dose},
    {"tip_time", &p.tip_time}, {"t_end", &p.t_end}, {"h", &p.h}, {"dt_out", &p.dt_out}};
}

static std::string trim(std::string s) {
  s.erase(0, s.find_first_not_of(" \t\r")); s.erase(s.find_last_not_of(" \t\r") + 1); return s;
}
static std::string errs_line(int ln, const std::string& m) { return "ERROR: line " + std::to_string(ln) + ": " + m; }

// Parse a stream of "name = value" lines into Params; every problem is appended to err with its
// line number. Used for both files (CLI) and in-memory text (Android/JNI), so behaviour is identical.
static void load_stream(std::istream& in, Params& p, std::vector<std::string>& err) {
  auto reg = registry(p);
  std::string line; int ln = 0;
  while (std::getline(in, line)) {
    ++ln;
    if (auto c = line.find('#'); c != std::string::npos) line.erase(c);
    auto eq = line.find('=');
    if (eq == std::string::npos) {
      if (!trim(line).empty()) err.push_back(errs_line(ln, "expected 'name = value'"));
      continue;
    }
    const std::string key = trim(line.substr(0, eq)), val = trim(line.substr(eq + 1));
    char* end = nullptr;
    const double v = std::strtod(val.c_str(), &end);
    if (val.empty() || *end != '\0' || !std::isfinite(v)) { err.push_back(errs_line(ln, "bad number for '" + key + "'")); continue; }
    bool found = false;
    for (auto& e : reg) if (key == e.name) { *e.ptr = v; found = true; }
    if (!found) err.push_back(errs_line(ln, "unknown parameter '" + key + "'"));
  }
}

// File wrapper around load_stream (the CLI uses it; the C API does not, hence maybe_unused).
[[maybe_unused]] static void load(const std::string& path, Params& p, std::vector<std::string>& err) {
  std::ifstream in(path);
  if (!in) { err.push_back("ERROR: cannot open " + path); return; }
  load_stream(in, p, err);
}

// Derived quantities. Infection-free equilibrium: T* = lambda/dT; C* = C_rest (chosen by the user);
// H* from Eq 10 at C*; and the basal cytokine synthesis s0 that makes C* an exact equilibrium of
// Eq 8 (solve dC/dt = 0 with no infection). So "C_rest" is an input and s0 follows from it.
static void derive(Params& p) {
  p.Ts = p.lambda / p.dT;
  p.Cs = p.c_rest;
  p.Hs = p.H0 + p.dH * hill(p.Cs, p.KH, (int)p.n_hill);
  p.s0 = p.dC * p.Cs / vagal(p, p.Hs) - p.gamma * hill(p.Cs, p.KC, 2);
}

// Reject physically meaningless input (errors) and flag risky-but-legal input (warnings).
// Warning example: if Eq 11 fails, the cytokine subsystem may be bistable (the "storm class").
// Order matters: derive() runs only if all basic checks pass, because it divides by parameters.
static void validate(Params& p, std::vector<std::string>& err) {
  auto pos = [&](const char* n, double v) { if (!(v > 0)) err.push_back(std::string("ERROR: ") + n + " must be > 0"); };
  pos("lambda", p.lambda); pos("d_T", p.dT); pos("beta_W", p.bW); pos("k_W", p.kW);
  pos("k_Co", p.kCo); pos("delta_W", p.dW); pos("delta_Co", p.dCo); pos("p_W", p.pW);
  pos("p_T", p.pT); pos("c_W", p.cW); pos("c_T", p.cT); pos("K_F", p.KF); pos("K_C", p.KC);
  pos("K_H", p.KH); pos("d_C", p.dC); pos("tau_H", p.tau_H); pos("C_rest", p.c_rest);
  pos("C_crit", p.C_crit); pos("H_crit", p.H_crit); pos("t_end", p.t_end); pos("h", p.h);
  pos("dt_out", p.dt_out);
  if (p.bT < 0 || p.dN < 0 || p.gamma < 0 || p.rhoI < 0 || p.rhoV < 0 || p.VW0 < 0 || p.tip_dose < 0)
    err.push_back("ERROR: beta_T, delta_N, gamma, rho_I, rho_V, V_W0, tip_dose must be >= 0");
  if (p.alpha < 0 || p.alpha > 1) err.push_back("ERROR: alpha must be in [0,1]");
  if (p.sigma_v < 0 || p.sigma_v >= 1) err.push_back("ERROR: sigma_v must be in [0,1)");
  if (p.eps_ab < 0 || p.eps_ab >= 1) err.push_back("ERROR: eps_Ab must be in [0,1)");
  if (p.eps_clr < 0 || p.eps_clr > 1) err.push_back("ERROR: eps_Clr must be in [0,1]");
  if (p.n_hill < 1 || p.n_hill > 4 || p.n_hill != std::floor(p.n_hill)) err.push_back("ERROR: n_hill must be an integer 1..4");
  if (p.ext_floor != 0 && p.ext_floor != 1) err.push_back("ERROR: ext_floor must be 0 or 1");
  if (p.Hmax <= p.H0) err.push_back("ERROR: H_max must exceed H_0");
  if (p.tip_time < 0 || p.tip_time > p.t_end) err.push_back("ERROR: tip_time must lie in [0, t_end]");
  if (p.dt_out < p.h) err.push_back("ERROR: dt_out must be >= h");
  if (!err.empty()) return;
  derive(p);
  if (p.s0 < 0) err.push_back("ERROR: C_rest too low for these parameters (basal synthesis s_0 would be negative)");
  if (3.0 * std::sqrt(3.0) / 8.0 * p.gamma / p.KC >= p.dC)
    err.push_back("WARNING: Eq 11 violated - cytokine subsystem may be bistable (storm-class regime)");
}

// Basic reproduction number of WT alone, Eq 12: R0 = bW*pW*T*/(cW*delta_eff), delta_eff at C*.
static double R0(const Params& p) {
  return p.bW * p.pW * p.Ts / (p.cW * (p.dW + p.dN * p.Cs / (p.KF + p.Cs)));  // Eq 12
}

// ---- One simulation: integrate from the infection-free state plus a WT inoculum, apply the TIP dose, record metrics.
// ------------------------------------------------------------------ one run
// A recorded trajectory: output times t[k] and the full state y[k] at those times.
struct Traj { std::vector<double> t; std::vector<State> y; };
// Summary metrics of one run: peak WT load and its time, peak cytokine/pulse, extinction flag and
// time, final WT, the four objective integrals, and ok=false if an invariant (finite, non-negative) broke.
struct Result {
  double peakV = 0, tPeakV = 0, peakC = 0, peakH = 0, tExt = -1, Vfinal = 0;
  std::array<double, 4> I{};  // JV, JC, JH, U
  bool extinct = false, ok = true;
  int reduced = 0;
};
// Objective J of Eq 14 = w1*JV + w2*JC + w3*JH + w4*U (smaller is better).
static double Jw(const Params& p, const Result& r) { return p.w1 * r.I[0] + p.w2 * r.I[1] + p.w3 * r.I[2] + p.w4 * r.I[3]; }

// Simulate one host. Steps: (1) infection-free start with V_W0 inoculum; (2) at t = tip_time add
// tip_dose to the free TIP pool; (3) march with RK4 between output times, splitting a segment so
// the dose lands exactly on tip_time; (4) after every step check invariants, apply the optional
// extinction floor, and update peaks. EXTINCTION means V_W < V_ext and stays below until t_end.
static Result run(const Params& p, Traj* tr, std::vector<std::string>& err) {
  Result r;
  State y{};
  y[iT] = p.Ts; y[iVW] = p.VW0; y[iC] = p.Cs; y[iH] = p.Hs;
  double t = 0, lastAbove = -1;
  bool dosed = false;
  // One-shot TIP dose: an instantaneous jump of V_T (e.g. an injection), applied once at tip_time.
  auto dose_if_due = [&] { if (!dosed && t >= p.tip_time - 1e-12) { y[iVT] += p.tip_dose; dosed = true; } };
  // Update running maxima and remember the last moment V_W was still at/above V_ext (for t_ext).
  auto track = [&] {
    if (y[iVW] > r.peakV) { r.peakV = y[iVW]; r.tPeakV = t; }
    r.peakC = std::max(r.peakC, y[iC]); r.peakH = std::max(r.peakH, y[iH]);
    if (y[iVW] >= p.v_ext) lastAbove = t;
  };
  // Integrate to time "to". The step is the smaller of h and the stability limit 2.5/rho, then
  // divided evenly so the segment ends exactly at "to". Invariants are checked every step.
  auto advance = [&](double to) {
    const double rho = jac_bound(p, y);
    const double hmax = std::min(p.h, 2.5 / std::max(rho, 1e-12));  // pre-flight stability bound
    const int n = std::max(1, (int)std::ceil((to - t) / hmax - 1e-9));
    if (hmax < p.h) ++r.reduced;
    const double hs = (to - t) / n;
    for (int s = 0; s < n; ++s) {
      rk4(p, y, hs);
      t += hs;
      for (int i = 0; i < N; ++i)
        if (!std::isfinite(y[i]) || y[i] < -1e-6 * std::max(1.0, std::fabs(p.Ts))) {
          err.push_back("ERROR: invariant broken (non-finite or negative state) at t = " + std::to_string(t));
          r.ok = false; return;
        }
      // Extinction floor: a continuous ODE never reaches exactly zero, so when fewer than one virion
      // and fewer than one infected cell per host remain, WT is set to zero (proxy for stochastic extinction).
      if (p.ext_floor > 0.5 && y[iVW] < p.v_ext &&
          y[iEW] + y[iIW] + y[iEC] + y[iIC] < p.v_ext) {  // WT-free: less than one virion/cell per host
        y[iVW] = y[iEW] = y[iIW] = y[iEC] = y[iIC] = 0;
      }
      track();
    }
    t = to;
  };
  // Save the current state when a trajectory is requested (screening passes nullptr to save memory).
  auto emit = [&] { if (tr) { tr->t.push_back(t); tr->y.push_back(y); } };
  dose_if_due(); track(); emit();
  const int nOut = (int)std::ceil(p.t_end / p.dt_out - 1e-9);
  for (int k = 1; k <= nOut && r.ok; ++k) {
    const double target = std::min(k * p.dt_out, p.t_end);
    if (!dosed && p.tip_time > t && p.tip_time < target) { advance(p.tip_time); dose_if_due(); }
    if (!r.ok) break;
    advance(target);
    dose_if_due();
    emit();
  }
  for (int i = 0; i < 4; ++i) r.I[i] = y[iJV + i];
  r.Vfinal = y[iVW];
  r.extinct = r.ok && y[iVW] < p.v_ext;
  if (r.extinct) r.tExt = std::max(lastAbove, 0.0);
  return r;
}


// ---- Plots written as plain SVG text (vector, no libraries, viewable in any browser). The reader
// of the future needs nothing but a text editor to regenerate or inspect the figures.
// ------------------------------------------------------------------ SVG plots (no dependencies)
// One curve on a plot: x/y points, colour, dashed flag and a legend name.
struct Series { std::vector<double> x, y; std::string color; bool dashed = false; std::string name; };
// Clean numeric formatting: fixed decimals for human-scale values, trailing zeros stripped.
static std::string fmt(double v) {
  char b[48];
  if (!std::isfinite(v)) { std::snprintf(b, sizeof b, "nan"); return b; }
  double a = std::fabs(v);
  if (a == 0.0) return "0";
  if (a >= 1e-3 && a < 1e4) {
    int dec = (a >= 100) ? 0 : (a >= 10) ? 1 : (a >= 1) ? 2 : 3;
    std::snprintf(b, sizeof b, "%.*f", dec, v);
    std::string s(b);
    if (s.find('.') != std::string::npos) {
      while (!s.empty() && s.back() == '0') s.pop_back();
      if (!s.empty() && s.back() == '.') s.pop_back();
    }
    return s.empty() ? "0" : s;
  }
  std::snprintf(b, sizeof b, "%.2e", v);
  return b;
}
static const double kNaN = std::numeric_limits<double>::quiet_NaN();

// ---- Figure style (one font family, fixed sizes for every figure) ---------------------------------
static const char* kFont = "Arial,Helvetica,sans-serif";
static const double kTitle = 21, kPanelTitle = 17, kLabel = 15, kTick = 13, kLegend = 14, kNote = 13;

// Minimal markup for real symbols:  *V*  = italic,  _{W}  = subscript,  ^{2}  = superscript.
// Returns the SVG <tspan> string and the plain text (used to estimate the text width).
struct Rich { std::string svg, plain; };
static Rich rich(const std::string& s, double fs) {
  Rich r; bool ital = false; double shift = 0, prev = 0; std::string seg;
  auto flush = [&](double size_pct, double sh, bool it) {
    if (seg.empty()) return;
    char buf[160];
    std::snprintf(buf, sizeof buf, "<tspan%s%s%s>", it ? " font-style='italic'" : "",
                  size_pct < 100 ? (" font-size='" + std::to_string((int)size_pct) + "%'").c_str() : "",
                  sh != prev ? (" dy='" + fmt(sh - prev) + "'").c_str() : "");
    r.svg += buf; r.svg += seg; r.svg += "</tspan>"; r.plain += seg; seg.clear(); prev = sh;
  };
  size_t i = 0; double pct = 100;
  while (i < s.size()) {
    char ch = s[i];
    if (ch == '*') { flush(pct, shift, ital); ital = !ital; ++i; continue; }
    if ((ch == '_' || ch == '^') && i + 1 < s.size() && s[i + 1] == '{') {
      flush(pct, shift, ital);
      const size_t e = s.find('}', i);
      shift = (ch == '_') ? 0.22 * fs : -0.38 * fs; pct = 75;
      seg = s.substr(i + 2, e - i - 2); flush(pct, shift, ital);
      shift = 0; pct = 100; i = e + 1; continue;
    }
    if (ch == '&') seg += "&amp;"; else if (ch == '<') seg += "&lt;"; else seg += ch;
    ++i;
  }
  flush(pct, shift, ital);
  if (prev != 0) r.svg += "<tspan dy='" + fmt(-prev) + "'></tspan>";
  return r;
}
// Conservative text-width estimate (DejaVu Sans is the widest common fallback for Arial).
static double tw(const std::string& plain, double fs) { return 0.60 * fs * (double)plain.size(); }

// ---- Layout checker: every text box is recorded; overlaps, canvas overflow and text lying on a data line are reported.
struct Layout {
  struct Box { double x0, y0, x1, y1; std::string what; bool note; };
  struct Seg { double x0, y0, x1, y1; };
  std::vector<Box> t; std::vector<Seg> data; double W = 0, H = 0;
  void box(double x0, double y0, double x1, double y1, const std::string& what, bool note = false) { t.push_back({x0, y0, x1, y1, what, note}); }
  void text(double x, double y, const char* anchor, const std::string& plain, double fs, bool note = false) {
    const double w = tw(plain, fs), x0 = !std::strcmp(anchor, "end") ? x - w : !std::strcmp(anchor, "middle") ? x - w / 2 : x;
    box(x0, y - 0.80 * fs, x0 + w, y + 0.25 * fs, plain, note);
  }
  int check(const std::string& name) const {
    int bad = 0;
    for (size_t i = 0; i < t.size(); ++i) {
      const Box& a = t[i];
      if (a.x0 < 0 || a.y0 < 0 || a.x1 > W || a.y1 > H) { ++bad; std::fprintf(stderr, "LAYOUT %s: '%s' leaves the canvas\n", name.c_str(), a.what.c_str()); }
      for (size_t j = i + 1; j < t.size(); ++j) {
        const Box& b = t[j];
        if (a.x0 < b.x1 && b.x0 < a.x1 && a.y0 < b.y1 && b.y0 < a.y1) { ++bad; std::fprintf(stderr, "LAYOUT %s: text '%s' overlaps '%s'\n", name.c_str(), a.what.c_str(), b.what.c_str()); }
      }
      if (!a.note) continue;
      for (const Seg& s : data) {   // annotation/legend text must stay off the data lines (4 px clearance)
        const double len = std::hypot(s.x1 - s.x0, s.y1 - s.y0); const int n = std::max(1, (int)(len / 2));
        for (int k = 0; k <= n; ++k) {
          const double px = s.x0 + (s.x1 - s.x0) * k / n, py = s.y0 + (s.y1 - s.y0) * k / n;
          if (px > a.x0 - 4 && px < a.x1 + 4 && py > a.y0 - 4 && py < a.y1 + 4) { ++bad; std::fprintf(stderr, "LAYOUT %s: '%s' lies on a data line\n", name.c_str(), a.what.c_str()); k = n + 1; break; }
        }
      }
    }
    std::fprintf(stderr, "layout check [%s]: %d problem(s) in %zu text items\n", name.c_str(), bad, t.size());
    return bad;
  }
};

// Round-number tick step (1, 2 or 5 times a power of ten) giving about 'target' intervals.
static double nice_step(double lo, double hi, int target) {
  const double raw = std::max(hi - lo, 1e-300) / target, mag = std::pow(10.0, std::floor(std::log10(raw))), r = raw / mag;
  return (r < 1.5 ? 1 : r < 3 ? 2 : r < 7 ? 5 : 10) * mag;
}
// Tick label with as many decimals as the step needs (0.04, 0.02, 0.01), thousands separators, true minus sign.
static std::string fmt_tick(double v, double step) {
  const int dec = std::max(0, (int)std::ceil(-std::log10(step) - 1e-9));
  char b[64]; std::snprintf(b, sizeof b, "%.*f", dec, std::fabs(v) < step * 1e-9 ? 0.0 : v);
  std::string s(b), sign;
  if (!s.empty() && s[0] == '-') { sign = "\xE2\x88\x92"; s.erase(0, 1); }
  const size_t dot = s.find('.'); size_t ip = dot == std::string::npos ? s.size() : dot;
  for (int k = (int)ip - 3; k > 0; k -= 3) s.insert((size_t)k, ",");
  return sign + s;
}
static std::string pow10_label(int k, double fs, std::string* plain) {
  const std::string e = k < 0 ? "\xE2\x88\x92" + std::to_string(-k) : std::to_string(k);
  if (plain) *plain = "10" + e;
  return "10<tspan dy='" + fmt(-0.38 * fs) + "' font-size='75%'>" + e + "</tspan><tspan dy='" + fmt(0.38 * fs) + "'></tspan>";
}

// Draw one framed panel (axes, ticks, labels, optional horizontal safety line, curves) at (X,Y) of size WxH.
// logy plots log10 of the data. Returns SVG fragments; no state is kept (lay, when given, receives the text boxes).
// Style: bottom+left spines only, light grid, round-number ticks, no top/right spines, title above, reference-line label in the right margin.
static std::string panel(double X, double Y, double W, double H, const std::string& title, const std::string& ylab,
                         const std::vector<Series>& ss, bool logy, double hline, const std::string& hlab, double xlo_ = kNaN,
                         Layout* lay = nullptr) {
  const double L = 96, R = 70, T = 36, B = 56;
  double x0 = 1e300, x1 = -1e300, y0 = 1e300, y1 = -1e300;
  auto ty = [&](double v) { return logy ? std::log10(std::max(v, 1e-4)) : v; };
  for (auto& s : ss) for (size_t i = 0; i < s.x.size(); ++i) {
    x0 = std::min(x0, s.x[i]); x1 = std::max(x1, s.x[i]); y0 = std::min(y0, ty(s.y[i])); y1 = std::max(y1, ty(s.y[i])); }
  if (!std::isnan(hline)) { y0 = std::min(y0, ty(hline)); y1 = std::max(y1, ty(hline)); }
  if (!std::isnan(xlo_)) x0 = xlo_;
  if (x1 <= x0) x1 = x0 + 1;
  double ystep = 1;
  if (logy) { y0 = std::floor(y0); y1 = std::ceil(y1); if (y1 - y0 < 1) y1 = y0 + 1; ystep = (y1 - y0 > 6) ? 2 : 1; }
  else { y0 = std::min(0.0, y0); y1 = y1 * 1.06 + 1e-12; ystep = nice_step(y0, y1, 5); y1 = std::ceil(y1 / ystep - 1e-9) * ystep; y0 = std::floor(y0 / ystep) * ystep; }
  const double xstep = nice_step(x0, x1, 6);
  const double pw = W - L - R, ph = H - T - B;
  auto px = [&](double v) { return X + L + (v - x0) / (x1 - x0) * pw; };
  auto pyr = [&](double raw) { return Y + T + (1 - (raw - y0) / (y1 - y0)) * ph; };
  auto py = [&](double v) { return pyr(ty(v)); };
  std::ostringstream o;
  auto txt = [&](double x, double y, const char* anc, double fs, const std::string& markup, const char* extra, bool note) {
    const Rich r = rich(markup, fs);
    o << "<text x='" << x << "' y='" << y << "' font-size='" << fs << "' text-anchor='" << anc << "' " << extra << ">" << r.svg << "</text>\n";
    if (lay) lay->text(x, y, anc, r.plain, fs, note);
  };
  txt(X + L, Y + 22, "start", kPanelTitle, title, "font-weight='bold' fill='#1a1a1a'", false);
  o << "<rect x='" << X + L << "' y='" << Y + T << "' width='" << pw << "' height='" << ph << "' fill='#fafafa'/>\n";
  for (double fx = std::ceil(x0 / xstep - 1e-9) * xstep; fx <= x1 + xstep * 1e-9; fx += xstep) {   // x ticks
    const double xp = px(fx);
    o << "<line x1='" << xp << "' y1='" << Y + T << "' x2='" << xp << "' y2='" << Y + T + ph << "' stroke='#e4e4e4' stroke-width='1'/>\n";
    o << "<line x1='" << xp << "' y1='" << Y + T + ph << "' x2='" << xp << "' y2='" << Y + T + ph + 5 << "' stroke='#444' stroke-width='1.2'/>\n";
    txt(xp, Y + T + ph + 21, "middle", kTick, fmt_tick(fx, xstep), "fill='#333'", false);
  }
  for (double fy = std::ceil(y0 / ystep - 1e-9) * ystep; fy <= y1 + ystep * 1e-9; fy += ystep) {   // y ticks
    const double yp = pyr(fy);
    o << "<line x1='" << X + L << "' y1='" << yp << "' x2='" << X + L + pw << "' y2='" << yp << "' stroke='#e4e4e4' stroke-width='1'/>\n";
    o << "<line x1='" << X + L - 5 << "' y1='" << yp << "' x2='" << X + L << "' y2='" << yp << "' stroke='#444' stroke-width='1.2'/>\n";
    if (logy) {
      std::string pl; const std::string m = pow10_label((int)std::lround(fy), kTick, &pl);
      o << "<text x='" << X + L - 9 << "' y='" << yp + 4.5 << "' font-size='" << kTick << "' text-anchor='end' fill='#333'>" << m << "</text>\n";
      if (lay) lay->text(X + L - 9, yp + 4.5, "end", pl, kTick);
    } else txt(X + L - 9, yp + 4.5, "end", kTick, fmt_tick(fy, ystep), "fill='#333'", false);
  }
  o << "<line x1='" << X + L << "' y1='" << Y + T + ph << "' x2='" << X + L + pw << "' y2='" << Y + T + ph << "' stroke='#222' stroke-width='1.6'/>\n";   // bottom + left spines only
  o << "<line x1='" << X + L << "' y1='" << Y + T << "' x2='" << X + L << "' y2='" << Y + T + ph << "' stroke='#222' stroke-width='1.6'/>\n";
  txt(X + L + pw / 2, Y + H - 12, "middle", kLabel, "time (days)", "fill='#333'", false);
  { const double cx = X + 18, cy = Y + T + ph / 2; const Rich r = rich(ylab, kLabel);
    o << "<text x='" << cx << "' y='" << cy << "' font-size='" << kLabel << "' text-anchor='middle' fill='#333' transform='rotate(-90 " << cx << " " << cy << ")'>" << r.svg << "</text>\n";
    if (lay) lay->box(cx - 12, cy - tw(r.plain, kLabel) / 2, cx + 4, cy + tw(r.plain, kLabel) / 2, ylab); }
  if (!std::isnan(hline)) {   // safety / threshold line; its label sits in the right margin, never on the data
    o << "<line x1='" << X + L << "' y1='" << py(hline) << "' x2='" << X + L + pw << "' y2='" << py(hline) << "' stroke='#A32D2D' stroke-width='1.5' stroke-dasharray='5 4'/>\n";
    txt(X + L + pw + 8, py(hline) + 4.5, "start", kNote, hlab, "fill='#A32D2D'", true);
  }
  for (auto& s : ss) {   // data series
    o << "<polyline fill='none' stroke='" << s.color << "' stroke-width='2.4' stroke-linejoin='round'" << (s.dashed ? " stroke-dasharray='8 5'" : "") << " points='";
    for (size_t i = 0; i < s.x.size(); ++i) {
      if (i) o << " ";
      o << px(s.x[i]) << "," << py(s.y[i]);
      if (lay && i) lay->data.push_back({px(s.x[i - 1]), py(s.y[i - 1]), px(s.x[i]), py(s.y[i])});
    }
    o << "'/>\n";
  }
  return o.str();
}

// Draw a legend of coloured line samples (wraps to a new row at maxw). Returns the height used.
static double legend(std::ostringstream& o, double x, double y, const std::vector<std::pair<std::string, std::pair<std::string, bool>>>& items,
                     double maxw = 1e9, Layout* lay = nullptr) {
  double cx = x, cy = y;
  for (auto& it : items) {
    const Rich r = rich(it.first, kLegend);
    const double w = 30 + tw(r.plain, kLegend);
    if (cx > x && cx + w > x + maxw) { cx = x; cy += 26; }
    o << "<line x1='" << cx << "' y1='" << cy << "' x2='" << cx + 24 << "' y2='" << cy << "' stroke='" << it.second.first << "' stroke-width='2.4'" << (it.second.second ? " stroke-dasharray='8 5'" : "") << "/>\n";
    o << "<text x='" << cx + 30 << "' y='" << cy + 5 << "' font-size='" << kLegend << "' fill='#222'>" << r.svg << "</text>\n";
    if (lay) lay->text(cx + 30, cy + 5, "start", r.plain, kLegend);
    cx += w + 22;
  }
  return cy - y + 26;
}

// Extract column idx (optionally idx+idx2, e.g. total of two compartments) of a trajectory as a Series.
static Series col(const Traj& t, int idx, const std::string& c, bool dashed, const std::string& name, int idx2 = -1) {
  Series s; s.color = c; s.dashed = dashed; s.name = name; s.x = t.t;
  for (auto& v : t.y) s.y.push_back(v[idx] + (idx2 >= 0 ? v[idx2] : 0.0));
  return s;
}

// Write the SVG header (white page, one font family) and register the canvas size with the checker.
static void svg_open(std::ostringstream& o, double W, double H) {
  o << "<svg xmlns='http://www.w3.org/2000/svg' width='" << W << "' height='" << H << "' viewBox='0 0 " << W << " " << H << "' font-family='" << kFont << "' xml:space='preserve'>\n";
  o << "<rect width='100%' height='100%' fill='white'/>\n";
}

// 2x2 dashboard: free virus (log), infected cells, cytokine vs C_crit, pulse vs H_crit.
// "b" is the optional untreated trajectory drawn dashed for comparison. Returned as a string
// so the Android app can display it without touching the file system.
static std::string traj_svg_str(const std::string& title, const Params& p, const Traj& a, const Traj* b) {
  std::ostringstream o;
  const double W = 600, H = 340, gap = 10, top = 46, CW = 2 * W + 24, CH = top + 2 * H + gap + 64;
  Layout lay; lay.W = CW; lay.H = CH;
  svg_open(o, CW, CH);
  { const Rich r = rich(title, kTitle); o << "<text x='14' y='30' font-size='" << kTitle << "' font-weight='bold' fill='#111'>" << r.svg << "</text>\n"; lay.text(14, 30, "start", r.plain, kTitle); }
  const std::string blue = "#185FA5", coral = "#D85A30", grey = "#6B6B6B";
  const double X1 = 8, X2 = 16 + W, Y1 = top, Y2 = top + H + gap;
  { std::vector<Series> s{col(a, iVW, blue, false, "V_W")};
    s.push_back(col(a, iVT, coral, false, "V_T"));
    if (b) s.push_back(col(*b, iVW, grey, true, "V_W untreated"));
    o << panel(X1, Y1, W, H, "(a) Free virus (log scale)", "copies/mL", s, true, p.v_ext, "*V*_{ext}", kNaN, &lay); }
  { std::vector<Series> s{col(a, iIW, blue, false, "I_W"), col(a, iIC, coral, false, "I_Co")};
    if (b) s.push_back(col(*b, iIW, grey, true, "I_W untreated"));
    o << panel(X2, Y1, W, H, "(b) Productively infected cells", "cells/mL", s, false, kNaN, "", kNaN, &lay); }
  { std::vector<Series> s{col(a, iC, blue, false, "C")};
    if (b) s.push_back(col(*b, iC, grey, true, "C untreated"));
    o << panel(X1, Y2, W, H, "(c) Cytokine (IL-6-like)", "cytokine (pg/mL)", s, false, p.C_crit, "*C*_{crit}", kNaN, &lay); }
  { std::vector<Series> s{col(a, iH, blue, false, "H")};
    if (b) s.push_back(col(*b, iH, grey, true, "H untreated"));
    o << panel(X2, Y2, W, H, "(d) Pulse rate", "pulse (beats/min)", s, false, p.H_crit, "*H*_{crit}", kNaN, &lay); }
  std::vector<std::pair<std::string, std::pair<std::string, bool>>> items{{"with TIP", {blue, false}}, {"TIP virions *V*_{T} / *I*_{Co}", {coral, false}}};
  if (b) items.push_back({"untreated", {grey, true}});
  legend(o, 8 + 96, Y2 + H + 22, items, CW - 120, &lay);
  o << "</svg>\n";
  lay.check("trajectory");
  return o.str();
}
// Same figure written to a file (CLI use).
static void write_traj_svg(const std::string& path, const std::string& title, const Params& p, const Traj& a, const Traj* b) {
  std::ofstream(path) << traj_svg_str(title, p, a, b);
}

// ---- Random numbers and threading helpers (self-contained so results do not depend on the platform library).
// ------------------------------------------------------------------ RNG and helpers
// xoshiro256** generator seeded by SplitMix64: tiny, fast, high quality, and fully specified,
// so the Latin-hypercube design is identical on every compiler and machine for a given seed.
struct Rng {
  uint64_t s[4];
  explicit Rng(uint64_t seed) {
    for (auto& x : s) { seed += 0x9e3779b97f4a7c15ULL; uint64_t z = seed;
      z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL; z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL; x = z ^ (z >> 31); }
  }
  uint64_t next() {
    auto rotl = [](uint64_t x, int k) { return (x << k) | (x >> (64 - k)); };
    const uint64_t r = rotl(s[1] * 5, 7) * 9, t = s[1] << 17;
    s[2] ^= s[0]; s[3] ^= s[1]; s[1] ^= s[2]; s[0] ^= s[3]; s[2] ^= t; s[3] = rotl(s[3], 45);
    return r;
  }
  double u01() { return (next() >> 11) * (1.0 / 9007199254740992.0); }
};

// Run fn(i) for i in [0,n) on a pool of threads that grab the next index from an atomic counter.
// Each i writes only its own result slot, so the answer cannot depend on thread count or timing
// (checked: 1-thread and 4-thread screens give byte-identical CSV output).
template <class F> static void parallel_for(size_t n, int threads, F fn) {
  std::atomic<size_t> next{0};
  std::vector<std::thread> pool;  // std::thread (not jthread): works on every NDK toolchain
  for (int t = 0; t < threads; ++t)
    pool.emplace_back([&] { for (;;) { const size_t i = next.fetch_add(1); if (i >= n) break; fn(i); } });
  for (auto& th : pool) th.join();
}

// ---- TIP-profile screening: sample the design space, evaluate every profile with run(), rank,
// refine, find the Pareto front, and measure which parameters matter.
// ------------------------------------------------------------------ screening
// One search dimension: parameter name, range [lo,hi], log-scale flag, and its byte offset inside Params.
struct Dim { std::string name; double lo, hi; bool lg; size_t off; };
// Screening settings read from screen.txt: dimensions plus sample count, seed, refinement effort, weights, output prefix.
struct ScreenCfg {
  std::vector<Dim> dims;
  int samples = 2000, refine_top = 3, nm_iters = 80, threads = 0, top_k = 10;
  double seed = 12345, wfactor = 2, w_ext = 5;
  std::string prefix = "screen";
};

// Parse screen.txt. Control keys (seed, samples, ...) set ScreenCfg; any registry name followed by
// "min max [log]" becomes a search dimension. OFFSET TRICK: the dimension stores the member's byte
// offset in Params so evaluate() can write a sampled value into a copy of Params without a
// big switch. This is valid only because every swept member is a double inside one standard-layout struct.
static void load_screen_stream(std::istream& in, ScreenCfg& c, Params& base, std::vector<std::string>& err) {
  auto reg = registry(base);
  std::string line; int ln = 0;
  while (std::getline(in, line)) {
    ++ln;
    if (auto k = line.find('#'); k != std::string::npos) line.erase(k);
    auto eq = line.find('=');
    if (eq == std::string::npos) { if (!trim(line).empty()) err.push_back(errs_line(ln, "expected 'name = value'")); continue; }
    const std::string key = trim(line.substr(0, eq));
    std::istringstream rest(line.substr(eq + 1));
    std::vector<std::string> tok; std::string w;
    while (rest >> w) tok.push_back(w);
    auto num = [&](const std::string& s, double& v) { char* e = nullptr; v = std::strtod(s.c_str(), &e); return !s.empty() && *e == '\0' && std::isfinite(v); };
    double v = 0;
    if (key == "prefix") { if (tok.size() == 1) c.prefix = tok[0]; else err.push_back(errs_line(ln, "prefix needs one word")); continue; }
    const std::map<std::string, double*> ctl = {{"seed", &c.seed}, {"weight_factor", &c.wfactor}, {"w_ext", &c.w_ext}};
    const std::map<std::string, int*> ictl = {{"samples", &c.samples}, {"refine_top", &c.refine_top}, {"nm_iters", &c.nm_iters},
                                              {"threads", &c.threads}, {"top_k", &c.top_k}};
    if (ctl.count(key)) { if (tok.size() == 1 && num(tok[0], v)) *ctl.at(key) = v; else err.push_back(errs_line(ln, "bad value for " + key)); continue; }
    if (ictl.count(key)) { if (tok.size() == 1 && num(tok[0], v)) *ictl.at(key) = (int)v; else err.push_back(errs_line(ln, "bad value for " + key)); continue; }
    const Entry* hit = nullptr;
    for (auto& e : reg) if (key == e.name) hit = &e;
    if (!hit) { err.push_back(errs_line(ln, "unknown parameter '" + key + "'")); continue; }
    double lo, hi;
    if (tok.size() < 2 || !num(tok[0], lo) || !num(tok[1], hi) || lo >= hi) { err.push_back(errs_line(ln, "need 'min max [log]' with min < max for " + key)); continue; }
    const bool lg = tok.size() > 2 && tok[2] == "log";
    if (lg && lo <= 0) { err.push_back(errs_line(ln, "log range needs min > 0 for " + key)); continue; }
    c.dims.push_back({key, lo, hi, lg, (size_t)((char*)hit->ptr - (char*)&base)});
  }
  if (c.dims.empty()) err.push_back("ERROR: screen file defines no search dimensions");
  if (c.samples < 10) err.push_back("ERROR: samples must be >= 10");
}

// File wrapper for load_screen_stream.
[[maybe_unused]] static void load_screen(const std::string& path, ScreenCfg& c, Params& base, std::vector<std::string>& err) {
  std::ifstream in(path);
  if (!in) { err.push_back("ERROR: cannot open " + path); return; }
  load_screen_stream(in, c, base, err);
}

// Result of evaluating one candidate profile: its parameter values x, the run Result, validity,
// feasibility (safety bounds respected), penalty, scalar objective F, and a class:
// 0 = feasible AND WT extinction (best), 1 = feasible, 2 = violates a safety bound.
struct Eval {
  std::vector<double> x;
  Result r;
  bool valid = false, feasible = false;
  double pen = 0, F = 0;
  int cls = 2;  // 0 = feasible + extinction, 1 = feasible, 2 = infeasible
};

// Map a unit-interval coordinate u in [0,1] to the real parameter range (linear or logarithmic).
static double dimval(const Dim& d, double u) {
  u = std::clamp(u, 0.0, 1.0);
  return d.lg ? std::exp(std::log(d.lo) + u * (std::log(d.hi) - std::log(d.lo))) : d.lo + u * (d.hi - d.lo);
}

// Evaluate one profile. F = J (Eq 14) + penalty, where the penalty is 1000x the relative excess of
// peak cytokine / pulse over their bounds plus w_ext*log10(final V_W / V_ext) so that "closer to
// extinction" is rewarded smoothly. Screening uses a coarser OUTPUT grid (dt_out >= 1 d) for
// speed; peaks are still tracked on every integration step, so results are not degraded.
static Eval evaluate(const Params& base, const ScreenCfg& c, const std::vector<double>& u) {
  Eval e;
  Params q = base;
  for (size_t j = 0; j < c.dims.size(); ++j) {
    const double v = dimval(c.dims[j], u[j]);
    *(double*)((char*)&q + c.dims[j].off) = v;
    e.x.push_back(v);
  }
  q.dt_out = std::max(q.dt_out, 1.0);  // coarser output grid in screening (peaks are tracked every step)
  std::vector<std::string> tmp;
  validate(q, tmp);
  for (auto& s : tmp) if (s.rfind("ERROR", 0) == 0) { e.F = 1e12; e.pen = 1e12; return e; }
  e.r = run(q, nullptr, tmp);
  if (!e.r.ok) { e.F = 1e12; e.pen = 1e12; return e; }
  e.valid = true;
  e.feasible = e.r.peakC <= q.C_crit && e.r.peakH <= q.H_crit;
  e.pen = 1e3 * std::max(0.0, e.r.peakC / q.C_crit - 1) + 1e3 * std::max(0.0, e.r.peakH / q.H_crit - 1)
        + c.w_ext * std::log10(std::max(e.r.Vfinal, q.v_ext) / q.v_ext);
  e.F = Jw(q, e.r) + e.pen;
  e.cls = !e.feasible ? 2 : (e.r.extinct ? 0 : 1);
  return e;
}

// Nelder-Mead downhill simplex (Nelder & Mead, 1965) in the unit cube, used to polish the best
// sampled profiles. Derivative-free because the objective has kinks (extinction floor, bounds).
// Standard moves: reflect (-1), expand (-2), contract (+-0.5), shrink (0.5). Coordinates are
// clamped to [0,1], and the loop stops early when the simplex has collapsed (relative F spread 1e-7).
// Nelder-Mead in the unit cube
static Eval nelder_mead(const Params& base, const ScreenCfg& c, std::vector<double> u0, int& evals) {
  const size_t d = u0.size();
  auto f = [&](std::vector<double>& u) { for (auto& v : u) v = std::clamp(v, 0.0, 1.0); ++evals; return evaluate(base, c, u); };
  std::vector<std::vector<double>> V(d + 1, u0);
  std::vector<Eval> E(d + 1);
  for (size_t i = 0; i < d; ++i) { V[i + 1][i] += (V[i + 1][i] > 0.9 ? -0.1 : 0.1); }
  for (size_t i = 0; i <= d; ++i) E[i] = f(V[i]);
  for (int it = 0; it < c.nm_iters; ++it) {
    std::vector<size_t> o(d + 1); std::iota(o.begin(), o.end(), 0);
    std::sort(o.begin(), o.end(), [&](size_t a, size_t b) { return E[a].F < E[b].F; });
    auto V2 = V; auto E2 = E;
    for (size_t i = 0; i <= d; ++i) { V[i] = V2[o[i]]; E[i] = E2[o[i]]; }
    if (std::fabs(E[d].F - E[0].F) < 1e-7 * (std::fabs(E[0].F) + 1e-9)) break;
    std::vector<double> cen(d, 0.0);
    for (size_t i = 0; i < d; ++i) for (size_t j = 0; j < d; ++j) cen[j] += V[i][j] / d;
    auto mix = [&](double t) { std::vector<double> v(d); for (size_t j = 0; j < d; ++j) v[j] = cen[j] + t * (V[d][j] - cen[j]); return v; };
    auto xr = mix(-1.0); Eval er = f(xr);
    if (er.F < E[0].F) { auto xe = mix(-2.0); Eval ee = f(xe); if (ee.F < er.F) { V[d] = xe; E[d] = ee; } else { V[d] = xr; E[d] = er; } }
    else if (er.F < E[d - 1].F) { V[d] = xr; E[d] = er; }
    else {
      auto xc = er.F < E[d].F ? mix(-0.5) : mix(0.5); Eval ec = f(xc);
      if (ec.F < std::min(er.F, E[d].F)) { V[d] = xc; E[d] = ec; }
      else for (size_t i = 1; i <= d; ++i) { for (size_t j = 0; j < d; ++j) V[i][j] = V[0][j] + 0.5 * (V[i][j] - V[0][j]); E[i] = f(V[i]); }
    }
  }
  size_t b = 0;
  for (size_t i = 1; i <= d; ++i) if (E[i].F < E[b].F) b = i;
  return E[b];
}

// Average ranks (ties share the mean rank): the rank transform that turns Pearson into Spearman/PRCC.
static std::vector<double> ranks(const std::vector<double>& v) {
  std::vector<size_t> o(v.size()); std::iota(o.begin(), o.end(), 0);
  std::sort(o.begin(), o.end(), [&](size_t a, size_t b) { return v[a] < v[b]; });
  std::vector<double> r(v.size());
  for (size_t i = 0; i < o.size();) {
    size_t j = i; while (j + 1 < o.size() && v[o[j + 1]] == v[o[i]]) ++j;
    for (size_t k = i; k <= j; ++k) r[o[k]] = 0.5 * (i + j) + 1;
    i = j + 1;
  }
  return r;
}
// Gauss-Jordan matrix inverse with partial pivoting; returns false if (near-)singular. Small matrices only.
static bool invert(std::vector<std::vector<double>>& A) {
  const size_t n = A.size();
  std::vector<std::vector<double>> I(n, std::vector<double>(n, 0.0));
  for (size_t i = 0; i < n; ++i) I[i][i] = 1;
  for (size_t c = 0; c < n; ++c) {
    size_t p = c; for (size_t r = c + 1; r < n; ++r) if (std::fabs(A[r][c]) > std::fabs(A[p][c])) p = r;
    if (std::fabs(A[p][c]) < 1e-12) return false;
    std::swap(A[p], A[c]); std::swap(I[p], I[c]);
    const double d = A[c][c];
    for (size_t j = 0; j < n; ++j) { A[c][j] /= d; I[c][j] /= d; }
    for (size_t r = 0; r < n; ++r) if (r != c) { const double f = A[r][c];
      for (size_t j = 0; j < n; ++j) { A[r][j] -= f * A[c][j]; I[r][j] -= f * I[c][j]; } }
  }
  A = I; return true;
}
// PRCC (Marino et al., 2008): rank-transform inputs and output, build the correlation matrix, invert
// it, and read partial correlations from the inverse. |PRCC| near 1 means the input strongly drives
// the output with the other inputs held fixed; sign gives direction.
// partial rank correlation of every input with output y (Marino et al. 2008)
static std::vector<double> prcc(const std::vector<std::vector<double>>& X, const std::vector<double>& y) {
  const size_t d = X.size(), n = y.size();
  std::vector<std::vector<double>> R;
  for (auto& col : X) R.push_back(ranks(col));
  R.push_back(ranks(y));
  std::vector<std::vector<double>> C(d + 1, std::vector<double>(d + 1, 0.0));
  std::vector<double> mean(d + 1, 0.0), sd(d + 1, 0.0);
  for (size_t i = 0; i <= d; ++i) { for (double v : R[i]) mean[i] += v / n; for (double v : R[i]) sd[i] += (v - mean[i]) * (v - mean[i]); sd[i] = std::sqrt(sd[i]); }
  for (size_t i = 0; i <= d; ++i) for (size_t j = 0; j <= d; ++j) {
    double s = 0; for (size_t k = 0; k < n; ++k) s += (R[i][k] - mean[i]) * (R[j][k] - mean[j]);
    C[i][j] = (sd[i] > 0 && sd[j] > 0) ? s / (sd[i] * sd[j]) : (i == j ? 1.0 : 0.0);
  }
  std::vector<double> out(d, kNaN);
  if (!invert(C)) return out;
  for (size_t j = 0; j < d; ++j) out[j] = -C[j][d] / std::sqrt(C[j][j] * C[d][d]);
  return out;
}

// Format a double with 4 significant digits for reports.
static std::string fv(double v) { char b[32]; std::snprintf(b, sizeof b, "%.4g", v); return b; }

// The whole screening workflow. Deterministic for a given seed (and thread count does not matter).
// Writes <prefix>_report.txt, _samples.csv, _pareto.svg and _best.svg. Returns 0 on success.
static int screen(Params base, const ScreenCfg& c, std::vector<std::string>& err) {
  const auto t0 = std::chrono::steady_clock::now();
  const int nt = c.threads > 0 ? c.threads : std::max(1u, std::thread::hardware_concurrency());
  const size_t n = (size_t)c.samples, d = c.dims.size();
  // 1. LATIN HYPERCUBE SAMPLING (McKay et al., 1979): each dimension is cut into n equal slices and
  // every slice is used exactly once (random permutation + jitter), covering the space evenly.
  // Latin hypercube
  Rng rng((uint64_t)c.seed);
  std::vector<std::vector<double>> U(n, std::vector<double>(d));
  for (size_t j = 0; j < d; ++j) {
    std::vector<size_t> perm(n); std::iota(perm.begin(), perm.end(), 0);
    for (size_t i = n - 1; i > 0; --i) std::swap(perm[i], perm[rng.next() % (i + 1)]);
    for (size_t i = 0; i < n; ++i) U[i][j] = (perm[i] + rng.u01()) / n;
  }
  std::vector<Eval> E(n);
  parallel_for(n, nt, [&](size_t i) { E[i] = evaluate(base, c, U[i]); });

  // 2. Untreated reference run (tip_dose = 0) used for comparison in the report and figures.
  // untreated reference
  Params ub = base; ub.tip_dose = 0; std::vector<std::string> tmp;
  Traj trU; Result ru = run(ub, &trU, tmp);

  // 3. Ranking: lower class first (A feasible+extinct, B feasible, C unsafe), then lower F.
  // Ties break on sample index so the order is total and reproducible.
  // ranking (class, then scalar F)
  std::vector<size_t> ord; for (size_t i = 0; i < n; ++i) if (E[i].valid) ord.push_back(i);
  auto better = [&](size_t a, size_t b) { if (E[a].cls != E[b].cls) return E[a].cls < E[b].cls; return E[a].F < E[b].F; };
  std::sort(ord.begin(), ord.end(), [&](size_t a, size_t b) { return better(a, b) || (!better(b, a) && a < b); });
  size_t nValid = ord.size(), nFeas = 0, nExt = 0;
  for (size_t i : ord) { nFeas += E[i].feasible; nExt += E[i].cls == 0; }
  if (nValid == 0) { err.push_back("ERROR: no valid samples - check the search bounds"); return 1; }

  // 4. Pareto front over (peak V_W, peak C, TIP exposure U) among feasible samples: keep a profile
  // only if no kept profile is at least as good in all three and strictly better in one.
  // Pareto front over (peak V_W, peak C, exposure U) among feasible samples
  std::vector<size_t> feas; for (size_t i : ord) if (E[i].feasible) feas.push_back(i);
  std::sort(feas.begin(), feas.end(), [&](size_t a, size_t b) { return E[a].r.peakV < E[b].r.peakV || (E[a].r.peakV == E[b].r.peakV && a < b); });
  std::vector<size_t> front;
  for (size_t i : feas) {
    bool dom = false;
    for (size_t f : front) if (E[f].r.peakC <= E[i].r.peakC && E[f].r.I[3] <= E[i].r.I[3] &&
        (E[f].r.peakC < E[i].r.peakC || E[f].r.I[3] < E[i].r.I[3] || E[f].r.peakV < E[i].r.peakV)) { dom = true; break; }
    if (!dom) front.push_back(i);
  }
  std::vector<char> onFront(n, 0); for (size_t f : front) onFront[f] = 1;

  // 5. Polish the best refine_top sampled profiles with Nelder-Mead.
  // Nelder-Mead refinement of the best few
  struct Ref { size_t from; Eval e; int evals; };
  std::vector<Ref> refined;
  for (int k = 0; k < c.refine_top && k < (int)ord.size(); ++k) {
    int ev = 0; Eval best = nelder_mead(base, c, U[ord[k]], ev);
    refined.push_back({ord[k], best, ev});
  }
  std::sort(refined.begin(), refined.end(), [](const Ref& a, const Ref& b) { return a.e.F < b.e.F; });

  // 6. Weight sensitivity: scale each objective weight by x and 1/x (weight_factor) and check whether
  // the baseline winner stays first; reports how robust the ranking is to subjective weights.
  // weight sensitivity: re-rank under scaled objective weights
  const double w0[4] = {base.w1, base.w2, base.w3, base.w4};
  auto score = [&](size_t i, const double* w) { return w[0] * E[i].r.I[0] + w[1] * E[i].r.I[1] + w[2] * E[i].r.I[2] + w[3] * E[i].r.I[3] + E[i].pen; };
  const size_t win = ord[0];
  const char* wn[4] = {"w1 (virus)", "w2 (cytokine)", "w3 (pulse)", "w4 (exposure)"};
  std::ostringstream ws; int same = 0, tot = 0, worst = 1;
  for (int k = 0; k < 4; ++k) for (int s = 0; s < 2; ++s) {
    double w[4] = {w0[0], w0[1], w0[2], w0[3]}; w[k] *= (s ? 1.0 / c.wfactor : c.wfactor);
    size_t top = ord[0]; int rankWin = 1;
    for (size_t i : ord) { auto lt = [&](size_t a, size_t b) { return E[a].cls != E[b].cls ? E[a].cls < E[b].cls : score(a, w) < score(b, w); };
      if (lt(i, top)) top = i;
      if (i != win && lt(i, win)) ++rankWin; }
    ++tot; same += (top == win); worst = std::max(worst, rankWin);
    ws << "  " << wn[k] << (s ? " / " : " x ") << fv(c.wfactor) << " : top-1 " << (top == win ? "unchanged" : "changed (sample " + std::to_string(top) + ")")
       << ", baseline winner rank " << rankWin << "\n";
  }

  // 7. Global sensitivity: PRCC of every swept parameter against log10 peak V_W, peak C and TIP exposure.
  // PRCC
  std::vector<std::vector<double>> X(d); std::vector<double> yV, yC, yU;
  for (size_t i : ord) { for (size_t j = 0; j < d; ++j) X[j].push_back(E[i].x[j]);
    yV.push_back(std::log10(std::max(E[i].r.peakV, 1e-4))); yC.push_back(E[i].r.peakC); yU.push_back(E[i].r.I[3]); }
  auto pV = prcc(X, yV), pC = prcc(X, yC), pU = prcc(X, yU);

  // 8. Outputs: text report, CSV of every sample (for any later analysis), and two SVG figures.
  // ---- report
  const std::string rep = c.prefix + "_report.txt";
  std::ofstream o(rep);
  o << "# TIP-profile screening report\n";
  o << "samples = " << n << ", valid = " << nValid << ", feasible (C<=C_crit, H<=H_crit) = " << nFeas
    << ", feasible + WT extinction = " << nExt << ", Pareto front size = " << front.size() << "\n";
  o << "seed = " << (uint64_t)c.seed << ", t_end = " << base.t_end << " d, V_ext = " << base.v_ext << "\n";
  o << "R0_WT = " << fv(R0(base)) << "\n";
  o << "untreated reference: peak_V_W = " << fv(ru.peakV) << ", peak_C = " << fv(ru.peakC) << ", peak_H = " << fv(ru.peakH)
    << ", extinct = " << (ru.extinct ? "yes" : "no") << ", J = " << fv(Jw(base, ru)) << "\n\n";
  o << "class: A = feasible + extinction (high priority), B = feasible, C = violates safety bound\n";
  auto head = [&] { o << "rank class"; for (auto& dm : c.dims) o << " " << dm.name; o << " peakV peakC peakH t_ext U_AUC_VT F\n"; };
  auto rowf = [&](size_t id, const Eval& e) {
    o << "(" << id << ") " << "ABC"[e.cls];
    for (double v : e.x) o << " " << fv(v);
    o << " " << fv(e.r.peakV) << " " << fv(e.r.peakC) << " " << fv(e.r.peakH) << " " << (e.r.extinct ? fv(e.r.tExt) : "none")
      << " " << fv(e.r.I[3]) << " " << fv(e.F) << "\n";
  };
  o << "== Top " << c.top_k << " profiles (sampled) ==\n"; head();
  for (int k = 0; k < c.top_k && k < (int)ord.size(); ++k) { o << (k + 1) << " "; rowf(ord[k], E[ord[k]]); }
  o << "\n== Nelder-Mead refinement of the top " << refined.size() << " ==\n"; head();
  for (auto& r : refined) { o << "from " << r.from << " (" << r.evals << " evals) "; rowf(r.from, r.e); o << "  F before = " << fv(E[r.from].F) << "\n"; }
  o << "\n== Pareto front (peak V_W, peak C, exposure), sorted by peak V_W; first 25 of " << front.size() << " ==\n"; head();
  for (size_t k = 0; k < front.size() && k < 25; ++k) { o << (k + 1) << " "; rowf(front[k], E[front[k]]); }
  o << "\n== PRCC (|value| near 1 = strong influence on the output) ==\nparameter  log10_peakV  peakC  exposure_U\n";
  for (size_t j = 0; j < d; ++j) o << c.dims[j].name << "  " << fv(pV[j]) << "  " << fv(pC[j]) << "  " << fv(pU[j]) << "\n";
  o << "\n== Weight sensitivity (each weight scaled by x" << fv(c.wfactor) << " and 1/" << fv(c.wfactor) << ") ==\n" << ws.str();
  o << "  baseline winner stays top-1 in " << same << " of " << tot << " scenarios (worst rank " << worst << ")\n";
  o << "\nNote: ranking is conditional on the illustrative parameter set; all outputs are simulated, not predictions.\n";
  o.close();

  std::ofstream cs(c.prefix + "_samples.csv");
  cs << "id"; for (auto& dm : c.dims) cs << "," << dm.name; cs << ",peakV,peakC,peakH,extinct,t_ext,U,F,class,pareto\n";
  for (size_t i = 0; i < n; ++i) {
    cs << i; for (size_t j = 0; j < d; ++j) cs << "," << E[i].x[j];
    cs << "," << E[i].r.peakV << "," << E[i].r.peakC << "," << E[i].r.peakH << "," << E[i].r.extinct << "," << E[i].r.tExt
       << "," << E[i].r.I[3] << "," << E[i].F << "," << "ABC"[E[i].cls] << "," << (int)onFront[i] << "\n";
  }

  // Figure: peak WT load (log) vs peak cytokine for all samples; Pareto points dark blue; red x = untreated.
  // Pareto scatter (same style as the dashboard: round-number ticks, light grid, bottom+left spines only, legend below the plot)
  { std::ostringstream s;
    const double W = 860, H = 520, L = 96, R = 80, T = 62, B = 66;
    double xl = 1e300, xh = -1e300, yh = base.C_crit * 1.05;
    auto lx = [&](double v) { return std::log10(std::max(v, 1e-4)); };
    for (size_t i : ord) { xl = std::min(xl, lx(E[i].r.peakV)); xh = std::max(xh, lx(E[i].r.peakV)); yh = std::max(yh, E[i].r.peakC * 1.02); }
    xl = std::min(xl, lx(ru.peakV)); xh = std::max(xh, lx(ru.peakV)); yh = std::max(yh, ru.peakC * 1.02);
    xl = std::floor(xl); xh = std::ceil(xh); if (xh <= xl) xh = xl + 1;
    yh = std::min(yh, base.C_crit * 2.0);
    const double ystep = nice_step(0, yh, 5); yh = std::ceil(yh / ystep - 1e-9) * ystep;
    const double xs = (xh - xl > 6) ? 2 : 1, pw = W - L - R, ph = H - T - B;
    auto px = [&](double v) { return L + (lx(v) - xl) / (xh - xl) * pw; };
    auto py = [&](double v) { return T + (1 - std::min(v, yh) / yh) * ph; };
    // legend entries (marker colour, label); wrapped onto rows that fit the canvas width
    const char* lg[5][2] = {{"#0C447C", "Pareto front"}, {"#1D9E75", "feasible + extinction"}, {"#85B7EB", "feasible"}, {"#B4B2A9", "violates a safety bound"}, {"#A32D2D", "untreated (x)"}};
    int rows = 1; { double cx = L; for (auto& it : lg) { const double w = 22 + tw(it[1], kLegend) + 28; if (cx > L && cx + w > W - R) { ++rows; cx = L; } cx += w; } }
    const double CH = H + 20 + 26 * rows;
    Layout lay; lay.W = W; lay.H = CH;
    svg_open(s, W, CH);
    auto txt = [&](double x, double y, const char* anc, double fs, const std::string& markup, const char* extra) {
      const Rich r = rich(markup, fs);
      s << "<text x='" << x << "' y='" << y << "' font-size='" << fs << "' text-anchor='" << anc << "' " << extra << ">" << r.svg << "</text>\n";
      lay.text(x, y, anc, r.plain, fs);
    };
    txt(L, 32, "start", kTitle, "TIP screening: peak WT viral load vs peak cytokine", "font-weight='bold' fill='#111'");
    s << "<rect x='" << L << "' y='" << T << "' width='" << pw << "' height='" << ph << "' fill='#fafafa'/>\n";
    for (double fx = xl; fx <= xh + 1e-9; fx += xs) {
      const double xp = L + (fx - xl) / (xh - xl) * pw;
      s << "<line x1='" << xp << "' y1='" << T << "' x2='" << xp << "' y2='" << (H - B) << "' stroke='#e4e4e4' stroke-width='1'/>\n";
      s << "<line x1='" << xp << "' y1='" << (H - B) << "' x2='" << xp << "' y2='" << (H - B + 5) << "' stroke='#444' stroke-width='1.2'/>\n";
      std::string pl; const std::string m = pow10_label((int)std::lround(fx), kTick, &pl);
      s << "<text x='" << xp << "' y='" << (H - B + 22) << "' font-size='" << kTick << "' text-anchor='middle' fill='#333'>" << m << "</text>\n";
      lay.text(xp, H - B + 22, "middle", pl, kTick);
    }
    for (double fy = 0; fy <= yh + ystep * 1e-9; fy += ystep) {
      const double yp = T + (1 - fy / yh) * ph;
      s << "<line x1='" << L << "' y1='" << yp << "' x2='" << (W - R) << "' y2='" << yp << "' stroke='#e4e4e4' stroke-width='1'/>\n";
      s << "<line x1='" << (L - 5) << "' y1='" << yp << "' x2='" << L << "' y2='" << yp << "' stroke='#444' stroke-width='1.2'/>\n";
      txt(L - 9, yp + 4.5, "end", kTick, fmt_tick(fy, ystep), "fill='#333'");
    }
    s << "<line x1='" << L << "' y1='" << (H - B) << "' x2='" << (W - R) << "' y2='" << (H - B) << "' stroke='#222' stroke-width='1.6'/>\n";
    s << "<line x1='" << L << "' y1='" << T << "' x2='" << L << "' y2='" << (H - B) << "' stroke='#222' stroke-width='1.6'/>\n";
    txt(L + pw / 2, H - 14, "middle", kLabel, "peak *V*_{W} (copies/mL, log scale)", "fill='#333'");
    { const double cx = 22, cy = T + ph / 2; const Rich r = rich("peak cytokine *C* (pg/mL)", kLabel);
      s << "<text x='" << cx << "' y='" << cy << "' font-size='" << kLabel << "' text-anchor='middle' fill='#333' transform='rotate(-90 " << cx << " " << cy << ")'>" << r.svg << "</text>\n";
      lay.box(cx - 12, cy - tw(r.plain, kLabel) / 2, cx + 4, cy + tw(r.plain, kLabel) / 2, "y label"); }
    s << "<line x1='" << L << "' y1='" << py(base.C_crit) << "' x2='" << (W - R) << "' y2='" << py(base.C_crit) << "' stroke='#A32D2D' stroke-width='1.5' stroke-dasharray='5 4'/>\n";
    txt(W - R + 8, py(base.C_crit) + 4.5, "start", kNote, "*C*_{crit}", "fill='#A32D2D'");
    for (size_t i : ord) { if (E[i].cls == 2 || onFront[i]) continue;
      s << "<circle cx='" << px(E[i].r.peakV) << "' cy='" << py(E[i].r.peakC) << "' r='3' fill='" << (E[i].cls == 0 ? "#1D9E75" : "#85B7EB") << "' fill-opacity='0.65'/>\n"; }
    for (size_t i : ord) if (E[i].cls == 2) s << "<circle cx='" << px(E[i].r.peakV) << "' cy='" << py(E[i].r.peakC) << "' r='2.4' fill='#B4B2A9' fill-opacity='0.5'/>\n";
    for (size_t f : front) s << "<circle cx='" << px(E[f].r.peakV) << "' cy='" << py(E[f].r.peakC) << "' r='4.8' fill='#0C447C' stroke='white' stroke-width='1.2'/>\n";
    s << "<path d='M" << px(ru.peakV) - 7 << " " << py(ru.peakC) - 7 << "L" << px(ru.peakV) + 7 << " " << py(ru.peakC) + 7 << "M" << px(ru.peakV) - 7 << " " << py(ru.peakC) + 7 << "L" << px(ru.peakV) + 7 << " " << py(ru.peakC) - 7 << "' stroke='#A32D2D' stroke-width='2.8'/>\n";
    double lxp = L, lyp = H + 14;
    for (auto& it : lg) {
      const double w = 22 + tw(it[1], kLegend) + 28;
      if (lxp > L && lxp + w > W - R) { lxp = L; lyp += 26; }
      s << "<circle cx='" << lxp + 6 << "' cy='" << lyp << "' r='4.8' fill='" << it[0] << "'/>\n";
      txt(lxp + 18, lyp + 5, "start", kLegend, it[1], "fill='#222'");
      lxp += w;
    }
    s << "</svg>\n";
    lay.check("pareto");
    std::ofstream(c.prefix + "_pareto.svg") << s.str(); }

  // Figure: trajectory of the best profile (refined if it improved on the sampled winner) vs untreated.
  // best profile trajectory vs untreated
  { const Eval& bst = refined.empty() || refined[0].e.F >= E[win].F ? E[win] : refined[0].e;
    Params q = base;
    for (size_t j = 0; j < d; ++j) *(double*)((char*)&q + c.dims[j].off) = bst.x[j];
    validate(q, tmp);
    Traj tb; run(q, &tb, tmp);
    write_traj_svg(c.prefix + "_best.svg", "Best screened TIP profile vs untreated", q, tb, &trU); }

  const double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  std::printf("screened %zu profiles on %d thread(s) in %.1f s: valid %zu, feasible %zu, feasible+extinct %zu, Pareto %zu\n",
              n, nt, secs, nValid, nFeas, nExt, front.size());
  std::printf("best sampled F = %.4g  (class %c)  | refined best F = %.4g\n", E[win].F, "ABC"[E[win].cls], refined.empty() ? E[win].F : refined[0].e.F);
  std::printf("wrote %s, %s_samples.csv, %s_pareto.svg, %s_best.svg\n", rep.c_str(), c.prefix.c_str(), c.prefix.c_str(), c.prefix.c_str());
  return 0;
}

// ---- Self-test of the integrator: run the same problem with step h, h/2, h/4; for a 4th-order
// method the error ratio e1/e2 is ~2^4, so order = log2(e1/e2) should lie in [3.8, 4.2].
// ------------------------------------------------------------------ self-test
// Integrate with a fixed number of equal steps (multiplied by "mult"), dosing at tip_time; no stability shrinking, so errors are pure truncation errors.
[[maybe_unused]] static State fixed_run(const Params& p, int mult) {
  State y{}; y[iT] = p.Ts; y[iVW] = p.VW0; y[iC] = p.Cs; y[iH] = p.Hs;
  const double seg1 = p.tip_time, seg2 = p.t_end - p.tip_time;
  if (seg1 > 0) { const int n = (int)std::ceil(seg1 / p.h) * mult; for (int i = 0; i < n; ++i) rk4(p, y, seg1 / n); }
  y[iVT] += p.tip_dose;
  if (seg2 > 0) { const int n = (int)std::ceil(seg2 / p.h) * mult; for (int i = 0; i < n; ++i) rk4(p, y, seg2 / n); }
  return y;
}
// Max relative difference between two states over the dynamic block.
[[maybe_unused]] static double dist(const State& a, const State& b) {
  double m = 0;
  for (int i = 0; i < kDyn; ++i) m = std::max(m, std::fabs(a[i] - b[i]) / std::max(1.0, std::fabs(b[i])));
  return m;
}


// NOTE FOR FUTURE PORTERS: the C API below is the stable boundary. A Kotlin/JNI (or Python ctypes,
// WebAssembly, ...) wrapper needs only these four functions; the model code above never changes.
// ====================================================================== C API (for Kotlin / JNI)
// All text in, text out: no files needed for single runs, so it works inside an Android sandbox.
// Params text = the same "name = value" lines as input.txt (any subset; missing names use defaults).
// Every returned string must be released with ebov_free().
// Copy a std::string to malloc memory the caller releases with ebov_free() (C-compatible ownership).
static char* dup_str(const std::string& s) {
  char* p = (char*)std::malloc(s.size() + 1);
  if (p) std::memcpy(p, s.c_str(), s.size() + 1);
  return p;
}
// Helpers: does the message list contain a hard ERROR (warnings do not count)? / join messages with newlines.
static bool has_err(const std::vector<std::string>& e) { for (auto& s : e) if (s.rfind("ERROR", 0) == 0) return true; return false; }
static std::string join(const std::vector<std::string>& e) { std::string o; for (auto& s : e) o += s + "\n"; return o; }

// Unmangled names so any language can bind to them.
extern "C" {
const char* ebov_version() { return "ebov-tip-engine 1.0 (C++20)"; }
void ebov_free(char* p) { std::free(p); }

// format: 0 = summary (key=value lines + ERROR/WARNING lines), 1 = trajectory table, 2 = SVG plot
// Single run from text parameters. format 0 = key=value summary, 1 = trajectory table (space
// separated, header line first), 2 = SVG dashboard. Returns ERROR lines instead if input is invalid.
char* ebov_run(const char* params_text, int format) {
  Params p; std::vector<std::string> err;
  std::istringstream in(params_text ? params_text : "");
  load_stream(in, p, err);
  if (!has_err(err)) validate(p, err);
  if (has_err(err)) return dup_str(join(err));
  Traj tr; Result r = run(p, &tr, err);
  const bool haveU = p.tip_dose > 0;
  Traj trU; Result ru;
  if (haveU) { Params ub = p; ub.tip_dose = 0; std::vector<std::string> tmp; ru = run(ub, &trU, tmp); }
  if (format == 2) return dup_str(traj_svg_str(haveU ? "EBOV-TIP trajectory (with TIP vs untreated)" : "EBOV within-host trajectory (untreated)", p, tr, haveU ? &trU : nullptr));
  if (format == 1) {
    std::ostringstream o; o << "time T E_W I_W V_W E_Co I_Co V_T C H J P\n";
    for (size_t k = 0; k < tr.t.size(); ++k) { const State& y = tr.y[k];
      o << tr.t[k]; for (int i = 0; i < 9; ++i) o << ' ' << y[i];
      o << ' ' << p.w1 * y[iJV] + p.w2 * y[iJC] + p.w3 * y[iJH] + p.w4 * y[iU] << ' ' << y[iP] << "\n"; }
    return dup_str(o.str());
  }
  std::ostringstream o;
  o << "R0_WT=" << R0(p) << "\npeak_V_W=" << r.peakV << "\nt_peak_V_W=" << r.tPeakV << "\npeak_C=" << r.peakC
    << "\npeak_H=" << r.peakH << "\nC_within_bound=" << (r.peakC <= p.C_crit) << "\nH_within_bound=" << (r.peakH <= p.H_crit)
    << "\nextinct=" << r.extinct << "\nt_extinct=" << (r.extinct ? r.tExt : -1.0) << "\nexposure_AUC_V_T=" << r.I[3]
    << "\nobjective_J=" << Jw(p, r) << "\n";
  if (haveU) o << "untreated_peak_V_W=" << ru.peakV << "\nuntreated_peak_C=" << ru.peakC << "\nuntreated_peak_H=" << ru.peakH << "\n";
  if (R0(p) <= 1) err.push_back("WARNING: R_0 <= 1 for WT alone; not an acute-infection regime");
  o << join(err);
  return dup_str(o.str());
}

// (Output files need a writable prefix path; the report text is also returned so the UI can show it at once.)
// Screening: writes <out_prefix>_report.txt, _samples.csv, _pareto.svg, _best.svg and returns the report text
// (or ERROR lines). out_prefix example on Android: "/data/data/<pkg>/files/screen".
char* ebov_screen(const char* params_text, const char* screen_text, const char* out_prefix) {
  Params p; std::vector<std::string> err;
  std::istringstream in(params_text ? params_text : "");
  load_stream(in, p, err);
  if (!has_err(err)) validate(p, err);
  ScreenCfg c;
  std::istringstream sin(screen_text ? screen_text : "");
  load_screen_stream(sin, c, p, err);
  if (has_err(err)) return dup_str(join(err));
  if (out_prefix && *out_prefix) c.prefix = out_prefix;
  if (screen(p, c, err) != 0) return dup_str(join(err));
  std::ifstream rep(c.prefix + "_report.txt");
  std::ostringstream o; o << rep.rdbuf();
  return dup_str(join(err) + o.str());
}
}  // extern "C"

// main() is compiled out when the file is built as a library for the API test or JNI (-DEBOV_NO_MAIN).
#ifndef EBOV_NO_MAIN
// Command-line driver. Modes: default (single run), --screen, --selftest. Argument handling is
// positional with defaults; all input problems are written to errors.txt and nothing is run on error.
int main(int argc, char** argv) {
  const std::string mode = argc > 1 ? argv[1] : "";
  const bool self = mode == "--selftest", scr = mode == "--screen";
  const int a0 = (self || scr) ? 2 : 1;
  auto arg = [&](int k, const char* def) { return argc > a0 + k ? std::string(argv[a0 + k]) : std::string(def); };
  const std::string scfg = scr ? arg(0, "screen.txt") : "";
  const std::string fin = scr ? arg(1, "input.txt") : arg(0, "input.txt");
  const std::string fout = (self || scr) ? "" : arg(1, "output.txt");
  const std::string ferr = scr ? arg(2, "errors.txt") : (self ? "errors.txt" : arg(2, "errors.txt"));

  Params p;
  std::vector<std::string> err;
  load(fin, p, err);
  if (err.empty()) validate(p, err);
  auto has_error = [&] { for (auto& s : err) if (s.rfind("ERROR", 0) == 0) return true; return false; };
  auto flush_err = [&] { std::ofstream e(ferr); for (auto& s : err) e << s << "\n"; };
  if (has_error()) { flush_err(); std::fprintf(stderr, "Input rejected, see %s\n", ferr.c_str()); return 1; }

  if (self) {
    State a = fixed_run(p, 1), b = fixed_run(p, 2), c = fixed_run(p, 4);
    const double e1 = dist(a, b), e2 = dist(b, c), order = std::log2(e1 / e2);
    std::printf("RK4 order estimate: %.3f  (e1=%.3e, e2=%.3e)  %s\n", order, e1, e2, (order >= 3.8 && order <= 4.2) ? "PASS" : "CHECK");
    flush_err(); return 0;
  }
  if (scr) {
    ScreenCfg c;
    load_screen(scfg, c, p, err);
    if (has_error()) { flush_err(); std::fprintf(stderr, "Screen file rejected, see %s\n", ferr.c_str()); return 1; }
    const int rc = screen(p, c, err);
    flush_err(); return rc;
  }

  Traj tr; Result r = run(p, &tr, err);
  std::ofstream out(fout);
  out << "# time T E_W I_W V_W E_Co I_Co V_T C H J\n";
  for (size_t k = 0; k < tr.t.size(); ++k) {
    const State& y = tr.y[k];
    char buf[400];
    std::snprintf(buf, sizeof buf, "%.6g %.9g %.9g %.9g %.9g %.9g %.9g %.9g %.9g %.9g %.9g\n", tr.t[k], y[0], y[1], y[2], y[3], y[4], y[5], y[6], y[7], y[8],
                  p.w1 * y[iJV] + p.w2 * y[iJC] + p.w3 * y[iJH] + p.w4 * y[iU]);
    out << buf;
  }
  const double r0 = R0(p);
  if (r0 <= 1) err.push_back("WARNING: R_0 <= 1 for WT alone; not an acute-infection regime");
  if (r.reduced > 0) err.push_back("WARNING: step reduced for stability in " + std::to_string(r.reduced) + " output intervals");
  Traj trU; Result ru; const bool haveU = p.tip_dose > 0;
  if (haveU) { Params ub = p; ub.tip_dose = 0; std::vector<std::string> tmp; ru = run(ub, &trU, tmp); }
  out << "# ---- summary ----\n# R0_WT = " << r0 << "\n# R_eff_(entry+clearance only) = " << (1 - p.eps_ab) * r0 / (1 + p.kappa * p.eps_clr) << "\n";
  out << "# peak_V_W = " << r.peakV << " at t = " << r.tPeakV << " d" << (haveU ? "   (untreated: " + fv(ru.peakV) + ")" : "") << "\n";
  out << "# peak_C = " << r.peakC << (r.peakC > p.C_crit ? "  EXCEEDS C_crit" : "  within C_crit") << (haveU ? "   (untreated: " + fv(ru.peakC) + ")" : "") << "\n";
  out << "# peak_H = " << r.peakH << (r.peakH > p.H_crit ? "  EXCEEDS H_crit" : "  within H_crit") << (haveU ? "   (untreated: " + fv(ru.peakH) + ")" : "") << "\n";
  out << "# extinction (V_W < V_ext) = " << (r.extinct ? "YES" : "NO"); if (r.extinct) out << " , time = " << r.tExt << " d";
  out << "\n# objective_J = " << Jw(p, r) << "\n# exposure_AUC_V_T = " << r.I[3] << "\n";
  out.close();
  std::string svg = fout; if (auto dot = svg.rfind('.'); dot != std::string::npos) svg.erase(dot); svg += ".svg";
  write_traj_svg(svg, haveU ? "EBOV-TIP trajectory (with TIP vs untreated)" : "EBOV within-host trajectory (untreated)", p, tr, haveU ? &trU : nullptr);
  std::printf("R0=%.3g  peakV=%.3g (t=%.2f)  peakC=%.3g  peakH=%.3g  extinct=%s  J=%.4g  -> %s, %s\n",
              r0, r.peakV, r.tPeakV, r.peakC, r.peakH, r.extinct ? "yes" : "no", Jw(p, r), fout.c_str(), svg.c_str());
  flush_err();
  return r.ok ? 0 : 2;
}
#endif  // EBOV_NO_MAIN
