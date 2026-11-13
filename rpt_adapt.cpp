// =====================================================================================
// rpt_adapt.cpp - Adaptive-immunity and vaccination extension of the RPT-EBOV engine
// (C++20, standard library only).                                  Author: Reinford Nyachio
//
// PURPOSE
//   Adds three adaptive-immune states to the ten-state core model of code.cpp: the B-cell /
//   antibody-secreting pool B, the neutralising antibody concentration A and the cytotoxic
//   T-lymphocyte pool Z. Vaccination is represented ONLY through the initial conditions
//   (A(0), B(0), Z(0)); the equations are identical for naive and vaccinated hosts.
//   With adaptive = 0 the three new states are frozen and the programme reproduces the core
//   model exactly (this is how the extension is cross-validated against code.cpp).
//
// BUILD / RUN
//   g++ -std=c++20 -O2 rpt_adapt.cpp -o rpt_adapt
//   ./rpt_adapt name=value ... [traj=out.csv]    one run; prints a one-line summary
//   Any parameter below can be overridden as name=value on the command line.
//
// STATE (13 dynamic states):  T EW IW EC IC P VW VT C H | B A Z
// UNITS: cells/mL, copies/mL, days; antibody in ug/mL; cytokine in pg/mL; pulse in beats/min.
// =====================================================================================
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <map>
#include <string>
#include <vector>

// Index names for the state vector (core states first, then the adaptive block).
enum Idx { iT, iEW, iIW, iEC, iIC, iP, iVW, iVT, iC, iH, iB, iA, iZ, N };
using State = std::array<double, N>;

// All parameters. Core values are those of the reference parameter set of the main text;
// the adaptive block holds illustrative values stated in Table 15 of the document.
struct Params {
  // ---- core (see Section 5) ----
  double lambda = 5e4, dT = 0.07, bW = 5e-7, bT = 5e-7, kW = 6, kCo = 6;
  double dW = 2, dCo = 2, pW = 500, pT = 1500, alpha = 0.2, cW = 10, cT = 10;
  double dN = 0.5, KF = 50, kappa = 1;
  double c_rest = 5, dC = 5, gamma = 3, KC = 50, rhoI = 0.005, rhoV = 1e-6, rhoT = 0;
  double sigma_v = 0.5, H0 = 70, Hmax = 130, dH = 50, KH = 100, n_hill = 2, tau_H = 0.02;
  double eps_ab = 0, eps_clr = 0;      // constant (non-immune) blocking / clearance effects
  double primed = 1, v_ext = 2e-4, ext_floor = 1;
  double guard = 1;                    // 1 = stability pre-flight shrinks the step; 0 = fixed step (convergence studies)
  double VW0 = 100, tip_dose = 0, tip_time = 0, t_end = 30, h = 0.005, dt_out = 0.25;
  // ---- adaptive block ----
  double adaptive = 1;                 // 1 = adaptive states evolve; 0 = frozen (core model)
  double KA = 10;                      // antibody half-effect concentration (ug/mL)
  double eps_max = 0.9;                // maximal antibody entry-blocking fraction
  double xiT = 1;                      // antibody cross-recognition of the TIP (0 = TIP invisible)
  double add_tip_clr = 0;              // 1 = antibody clearance of TIP is ADDITIVE (worst case)
  double piB = 0.5, KB = 1e3, Bmax = 200, mu_decay = 0.005, B_naive = 1;
  double pA = 0.0165, dA = 0.033;      // pA is calibrated so that B = 100 sustains A = 50 ug/mL
  double piZ = 0.4, KZ = 1e3, Zmax = 1e4, muZd = 0.005, Z_naive = 1, muZ = 1e-5;
  double A0 = 0, B0 = -1, Z0 = -1;     // initial adaptive state (B0, Z0 < 0: derived from A0)
  // ---- optional initial perturbations (equilibrium test) ----
  double T0f = 1, C0 = -1, H0i = -1;
  // ---- derived ----
  double Ts = 0, Cs = 0, Hs = 0, s0 = 0;
};

static double hill(double x, double K, int n) {
  const double xn = std::pow(std::max(x, 0.0), n), Kn = std::pow(K, n);
  return xn / (Kn + xn);
}
// Vagal gain m(H), Eq 9.
static double vagal(const Params& p, double H) {
  return 1.0 - p.sigma_v * std::clamp((p.Hmax - H) / (p.Hmax - p.H0), 0.0, 1.0);
}
// Saturating antibody effect a(A) = A/(K_A + A) in [0,1): drives eps_Ab(A) and eps_Clr(A) (Eqs 17-18).
static double aeff(const Params& p, double A) { return std::max(A, 0.0) / (p.KA + std::max(A, 0.0)); }

// Right-hand side: Eqs 1-10 of the core model with antibody-dependent entry blocking and clearance,
// CTL killing, and the adaptive equations 19-21 (see the document, Section 11).
static void rhs(const Params& p, const State& y, State& d) {
  const double a = aeff(p, y[iA]) * p.adaptive;
  const double epsAb = 1.0 - (1.0 - p.eps_ab) * (1.0 - p.eps_max * a);   // combined entry blocking (Eq 18)
  const double epsClr = p.eps_clr + a * (1.0 - p.eps_clr);                // combined clearance effect (Eq 18)
  const double ent = 1.0 - epsAb;
  const double cyt = p.dN * y[iC] / (p.KF + y[iC]);
  const double ctl = p.adaptive * p.muZ * y[iZ];                           // CTL killing rate of infected cells
  const double inf = ent * p.bW * y[iT] * y[iVW];
  const double tip = ent * p.bT * y[iVT];
  const double supW = ent * p.bW * y[iVW];
  d[iT]  = p.lambda - p.dT * y[iT] - inf - tip * y[iT];
  d[iEW] = inf - p.kW * y[iEW] - tip * y[iEW];
  d[iIW] = p.kW * y[iEW] - (p.dW + cyt + ctl) * y[iIW] - tip * y[iIW];
  d[iEC] = tip * y[iEW] + p.primed * supW * y[iP] - p.kCo * y[iEC];
  d[iIC] = p.kCo * y[iEC] + tip * y[iIW] - (p.dCo + cyt + ctl) * y[iIC];
  d[iP]  = p.primed * (tip * y[iT] - supW * y[iP]) - p.dT * y[iP];
  d[iVW] = p.pW * y[iIW] + p.alpha * p.pW * y[iIC] - p.cW * (1.0 + p.kappa * epsClr) * y[iVW];
  // TIP clearance: multiplicative boost (default, same structure as WT) or additive worst case.
  const double cTeff = p.add_tip_clr > 0.5 ? p.cT + p.xiT * p.kappa * epsClr * p.cW
                                           : p.cT * (1.0 + p.xiT * p.kappa * epsClr);
  d[iVT] = p.pT * y[iIC] - cTeff * y[iVT];
  d[iC]  = vagal(p, y[iH]) * (p.s0 + p.rhoI * (y[iIW] + y[iIC]) + p.rhoV * y[iVW] + p.rhoT * y[iVT]
         + p.gamma * hill(y[iC], p.KC, 2)) - p.dC * y[iC];
  d[iH]  = (p.H0 + p.dH * hill(y[iC], p.KH, (int)p.n_hill) - y[iH]) / p.tau_H;
  // ---- adaptive block (Eqs 19-21) ----
  const double Itot = y[iIW] + y[iIC];
  d[iB] = p.adaptive * (p.piB * y[iB] * y[iVW] / (p.KB + y[iVW]) * (1.0 - y[iB] / p.Bmax)
                        - p.mu_decay * (y[iB] - p.B_naive));
  d[iA] = p.adaptive * (p.pA * y[iB] - p.dA * y[iA]);
  d[iZ] = p.adaptive * (p.piZ * y[iZ] * Itot / (p.KZ + Itot) * (1.0 - y[iZ] / p.Zmax)
                        - p.muZd * (y[iZ] - p.Z_naive));
}

// Classical fixed-step RK4.
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

// Infection-free equilibrium and the basal cytokine synthesis that makes C* exact (as in the core).
static void derive(Params& p) {
  p.Ts = p.lambda / p.dT;
  p.Cs = p.c_rest;
  p.Hs = p.H0 + p.dH * hill(p.Cs, p.KH, (int)p.n_hill);
  p.s0 = p.dC * p.Cs / vagal(p, p.Hs) - p.gamma * hill(p.Cs, p.KC, 2);
  if (p.B0 < 0) p.B0 = std::max(p.B_naive, p.A0 * p.dA / p.pA);   // memory pool consistent with A0
  if (p.Z0 < 0) p.Z0 = std::max(p.Z_naive, p.B0);                 // CTL memory scaled alike
}
// Basic reproduction number (Eq 12) and its antibody-modified form (Eq 13) at antibody level A.
static double R0(const Params& p) { return p.bW * p.pW * p.Ts / (p.cW * (p.dW + p.dN * p.Cs / (p.KF + p.Cs))); }
static double Reff(const Params& p, double A) {
  const double a = aeff(p, A);
  const double eAb = 1.0 - (1.0 - p.eps_ab) * (1.0 - p.eps_max * a), eCl = p.eps_clr + a * (1.0 - p.eps_clr);
  return (1.0 - eAb) * R0(p) / (1.0 + p.kappa * eCl);
}

// STABILITY GUARD. RK4 needs h*rho < ~2.78 (real axis), where rho is the spectral radius of the Jacobian J.
// J is formed by forward differences and rho is estimated by power iteration (geometric mean of the growth
// ratio over the last iterations, which also copes with complex eigenvalue pairs). A plain matrix-norm bound
// is far too pessimistic here (units differ by 10 orders of magnitude, e.g. production p_T = 2877 against a
// state of order 1), so the true spectral radius is estimated instead and a 25% safety factor is applied.
// advance() shrinks the step when h*rho would exceed 2.5; this matters for very large TIP doses, where
// beta_T*V_T reaches hundreds per day.
static double jac_bound(const Params& p, const State& y) {
  State f0, f1, z = y;
  rhs(p, y, f0);
  double J[N][N];
  for (int j = 0; j < N; ++j) {
    const double e = 1e-6 * std::max(1.0, std::fabs(y[j]));
    z[j] = y[j] + e;
    rhs(p, z, f1);
    z[j] = y[j];
    for (int i = 0; i < N; ++i) J[i][j] = (f1[i] - f0[i]) / e;
  }
  std::array<double, N> v;
  v.fill(1.0 / std::sqrt((double)N));
  double logsum = 0;
  int cnt = 0;
  for (int it = 0; it < 300; ++it) {
    std::array<double, N> w{};
    for (int i = 0; i < N; ++i) for (int j = 0; j < N; ++j) w[i] += J[i][j] * v[j];
    double nrm = 0;
    for (int i = 0; i < N; ++i) nrm += w[i] * w[i];
    nrm = std::sqrt(nrm);
    if (nrm < 1e-300) return 0.0;
    if (it >= 250) { logsum += std::log(nrm); ++cnt; }
    for (int i = 0; i < N; ++i) v[i] = w[i] / nrm;
  }
  return 1.25 * std::exp(logsum / cnt);
}

struct Result { double peakV = 0, tPeak = 0, peakC = 0, peakH = 0, Vfinal = 0, tExt = -1; bool extinct = false, ok = true; };

// One run. Writes the trajectory (one row per output time) to 'out' when it is non-null.
static Result run(const Params& p, std::FILE* out) {
  Result r;
  State y{};
  y[iT] = p.Ts * p.T0f; y[iVW] = p.VW0;
  y[iC] = p.C0 >= 0 ? p.C0 : p.Cs;
  y[iH] = p.H0i >= 0 ? p.H0i : p.Hs;
  y[iB] = p.adaptive > 0.5 ? p.B0 : p.B_naive; y[iA] = p.adaptive > 0.5 ? p.A0 : 0; y[iZ] = p.adaptive > 0.5 ? p.Z0 : p.Z_naive;
  double t = 0, lastAbove = -1;
  bool dosed = false;
  auto track = [&] {
    if (y[iVW] > r.peakV) { r.peakV = y[iVW]; r.tPeak = t; }
    r.peakC = std::max(r.peakC, y[iC]); r.peakH = std::max(r.peakH, y[iH]);
    if (y[iVW] >= p.v_ext) lastAbove = t;
  };
  auto dose = [&] { if (!dosed && t >= p.tip_time - 1e-12) { y[iVT] += p.tip_dose; dosed = true; } };
  auto emit = [&] {
    if (!out) return;
    std::fprintf(out, "%.6g", t);
    for (int i = 0; i < N; ++i) std::fprintf(out, ",%.15g", y[i]);
    std::fprintf(out, "\n");
  };
  auto advance = [&](double to) {
    const double rho = jac_bound(p, y);
    const double hmax = p.guard > 0.5 ? std::min(p.h, 2.5 / std::max(rho, 1e-12)) : p.h;   // pre-flight stability bound
    const int n = std::max(1, (int)std::ceil((to - t) / hmax - 1e-9));
    const double hs = (to - t) / n;
    for (int s = 0; s < n; ++s) {
      rk4(p, y, hs); t += hs;
      for (int i = 0; i < N; ++i)
        if (!std::isfinite(y[i]) || y[i] < -1e-6 * p.Ts) { r.ok = false; return; }
      if (p.ext_floor > 0.5 && y[iVW] < p.v_ext && y[iEW] + y[iIW] + y[iEC] + y[iIC] < p.v_ext)
        y[iVW] = y[iEW] = y[iIW] = y[iEC] = y[iIC] = 0;
      track();
    }
    t = to;
  };
  if (out) std::fprintf(out, "t,T,EW,IW,EC,IC,P,VW,VT,C,H,B,A,Z\n");
  dose(); track(); emit();
  const int nOut = (int)std::ceil(p.t_end / p.dt_out - 1e-9);
  for (int k = 1; k <= nOut && r.ok; ++k) {
    const double target = std::min(k * p.dt_out, p.t_end);
    if (!dosed && p.tip_time > t && p.tip_time < target) { advance(p.tip_time); dose(); }
    if (!r.ok) break;
    advance(target); dose(); emit();
  }
  r.Vfinal = y[iVW];
  r.extinct = r.ok && y[iVW] < p.v_ext;
  if (r.extinct) r.tExt = std::max(lastAbove, 0.0);
  return r;
}

// Command-line parameter table: name -> member.
static std::map<std::string, double*> table(Params& p) {
  return {{"lambda", &p.lambda}, {"d_T", &p.dT}, {"beta_W", &p.bW}, {"beta_T", &p.bT}, {"k_W", &p.kW}, {"k_Co", &p.kCo},
          {"delta_W", &p.dW}, {"delta_Co", &p.dCo}, {"p_W", &p.pW}, {"p_T", &p.pT}, {"alpha", &p.alpha}, {"c_W", &p.cW},
          {"c_T", &p.cT}, {"delta_N", &p.dN}, {"K_F", &p.KF}, {"kappa", &p.kappa}, {"C_rest", &p.c_rest}, {"d_C", &p.dC},
          {"gamma", &p.gamma}, {"K_C", &p.KC}, {"rho_I", &p.rhoI}, {"rho_V", &p.rhoV}, {"rho_T", &p.rhoT},
          {"sigma_v", &p.sigma_v}, {"H_0", &p.H0}, {"H_max", &p.Hmax}, {"Delta_H", &p.dH}, {"K_H", &p.KH}, {"n_hill", &p.n_hill},
          {"tau_H", &p.tau_H}, {"eps_Ab", &p.eps_ab}, {"eps_Clr", &p.eps_clr}, {"primed", &p.primed}, {"V_ext", &p.v_ext},
          {"ext_floor", &p.ext_floor}, {"guard", &p.guard}, {"V_W0", &p.VW0}, {"tip_dose", &p.tip_dose}, {"tip_time", &p.tip_time},
          {"t_end", &p.t_end}, {"h", &p.h}, {"dt_out", &p.dt_out}, {"adaptive", &p.adaptive}, {"K_A", &p.KA},
          {"eps_max", &p.eps_max}, {"xi_T", &p.xiT}, {"add_tip_clr", &p.add_tip_clr}, {"pi_B", &p.piB}, {"K_B", &p.KB},
          {"B_max", &p.Bmax}, {"mu_decay", &p.mu_decay}, {"B_naive", &p.B_naive}, {"p_A", &p.pA}, {"d_A", &p.dA},
          {"pi_Z", &p.piZ}, {"K_Z", &p.KZ}, {"Z_max", &p.Zmax}, {"mu_Zdecay", &p.muZd}, {"Z_naive", &p.Z_naive},
          {"mu_Z", &p.muZ}, {"A0", &p.A0}, {"B0", &p.B0}, {"Z0", &p.Z0}, {"T0_factor", &p.T0f}, {"C0", &p.C0}, {"H0_init", &p.H0i}};
}

int main(int argc, char** argv) {
  Params p;
  auto tab = table(p);
  std::string trajFile;
  for (int i = 1; i < argc; ++i) {
    const char* eq = std::strchr(argv[i], '=');
    if (!eq) { std::fprintf(stderr, "bad argument: %s\n", argv[i]); return 2; }
    const std::string k(argv[i], eq - argv[i]), v(eq + 1);
    if (k == "traj") { trajFile = v; continue; }
    auto it = tab.find(k);
    if (it == tab.end()) { std::fprintf(stderr, "unknown parameter: %s\n", k.c_str()); return 2; }
    *it->second = std::atof(v.c_str());
  }
  derive(p);
  std::FILE* out = trajFile.empty() ? nullptr : std::fopen(trajFile.c_str(), "w");
  const Result r = run(p, out);
  if (out) std::fclose(out);
  std::printf("ok=%d peakV=%.9g tPeak=%.5g peakC=%.9g peakH=%.9g Vfinal=%.9g extinct=%d tExt=%.5g R0=%.6g Reff0=%.6g\n",
              r.ok, r.peakV, r.tPeak, r.peakC, r.peakH, r.Vfinal, r.extinct, r.tExt, R0(p), Reff(p, p.A0));
  return r.ok ? 0 : 1;
}
