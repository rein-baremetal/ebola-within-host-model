# RPT-EBOV Advanced Theoretical Framework – Equation Specification

Independent Analytical Specification | Multi-Stage Entry, Host Interventions, Condensate Dynamics & Latency  
Aligned with core RPT-EBOV document (Nyachio, November 2026) and numerical results therein.

All equations below follow the Advanced Theoretical Framework PDF, with notation matching the core engine (code_annotated.cpp). Numerical illustrations (critical entry fractions, IC50 multiples) are retained only where they are consistent with the core document’s R₀ ≈ 8.7–8.9 and the analytic entry-blocking threshold 1 − 1/R₀ ≈ 0.885–0.888.

---

## 1. Overview & Staged Viral Entry Sub-Model

The core model uses a single-step infection term  
β_W · (1 − ε_Ab) · T · V_W.  

This extension expands it into a multi-step endosomal pipeline (proteolytic cleavage → NPC1 engagement → fusion priming):

```
dB/dt  = β_eff · V_W · T  − (k_c + d_B) · B
dK/dt  = k_c · B          − (k_f + d_B) · K
dE_W/dt = k_f · K         − k_W · E_W
```

- B  : bound endosomal cells  
- K  : protease-cleaved, NPC1-engaged states  
- k_c : cathepsin cleavage rate (day⁻¹)  
- k_f = φ_variant · k_f0 : fusion-priming rate  
- φ_variant ≥ 1 : glycoprotein variant entry-fusion multiplier (e.g. A82V); baseline = 1  

Effective entry rate:

```
β_eff = β_W · (1 − ε_Ab) · (1 − ε_drug(D))
```

Accelerated-fusion variants (φ_variant > 1) increase the basic reproductive ratio R₀ relative to the core reference and therefore elevate the critical entry-blocking fraction required for extinction.  

Core-document reference: pure entry blockade must exceed ≈ 88.5–88.8 % (1 − 1/R₀).  
Illustration retained from framework: a 4-fold acceleration can push the required blockade toward ≈ 97 % under otherwise identical parameters. This is a structural prediction, not a calibrated claim.

---

## 2. Host-Targeted Pharmacological Terms

Small-molecule entry inhibitors acting on cellular trafficking (e.g. TPC2 channels) affect wild-type (V_W) and therapeutic interfering particles (V_T) symmetrically through a Hill response:

```
ε_drug(D) = D^h / (IC₅₀^h + D^h)
```

- IC₅₀ : half-maximal concentration (framework reference value 51.7 μM for laudanosine-class compounds is illustrative)  
- h ∈ [0.5, 2.0] : Hill coefficient  

To reach the core-document baseline of ≈ 89 % entry reduction:

- h = 1.0 → D ≈ 8.1 × IC₅₀  
- h = 2.0 → D ≈ 2.8 × IC₅₀  

These multiples are pure Hill-function arithmetic and are independent of the viral dynamics; they remain valid under the core model’s entry-blocking formulation.

---

## 3. Condensate & Primed-Cell Dynamics

Viral condensates (NP/VP35 inclusion bodies) partition host mRNAs and viral genomes, altering interferon-stimulated gene (ISG) signalling.

Primed-cell superinfection (already present in the core model as the P compartment) is modulated by an interferon-protection factor π ∈ [0, 1]:

```
dE_Co/dt = β_W · (1 − π) · P · V_W − k_Co · E_Co
```

(π = 0 recovers the core primed-cell equation exactly.)

Cytokine production receives an additional condensate-driven gain:

```
… + g_cond · (I_W + I_Co) · ISG_factor(π)
```

where g_cond is a dimensionless coupling factor and ISG_factor(π) is a decreasing function of protection (exact functional form left open for calibration).

---

## 4. Counter-Regulatory Cytokine Subsystem & Latent Reservoirs

To prevent unrestrained inflammatory escalation a counter-regulatory mediator R (tissue-protective axes such as IL-10/IL-22) balances the pro-inflammatory cytokine C:

```
dC/dt = [ρ_I·(I_W + I_Co) + ρ_V·V_W + ρ_T·V_T + g_cond·I_total] / (1 + R/K_R) − d_C · (C − C₀)

dR/dt = σ_R · C / (K_σ + C) − d_R · R
```

A non-productive latent compartment L accounts for quiescent reservoirs:

```
dL/dt = v · (E_W + I_W) − (η + d_L) · L
```

- v  : latency seeding rate  
- η  : reactivation rate  
- d_L : latent-cell death  

**Extinction definitions (refined)**

- Productive extinction : active viraemia V_W < V_ext (the quantity already tracked by the core engine).  
- Sterilising extinction : both V_W < V_ext and latent pool L < L_ext.  

Because TIPs require active co-infection to replicate, cells in L remain inaccessible during latency. Consequently sterilising extinction may require structured multi-dose TIP schedules—an extension of the core persistence requirement (c_T ≪ c_W).

---

## 5. Model-Extension Parameter Reference

| Symbol      | Description                                      | Nominal / Range                          |
|-------------|--------------------------------------------------|------------------------------------------|
| φ_variant   | Glycoprotein variant entry-fusion multiplier     | ≥ 1.0 (baseline = 1.0)                   |
| k_c, k_f    | Cathepsin cleavage & NPC1-fusion priming rates   | Fitted constants (day⁻¹)                 |
| IC₅₀, h     | Half-maximal drug conc. & Hill coefficient       | illustrative 51.7 μM; h ∈ [0.5, 2.0]     |
| π           | Interferon protection factor on TIP-primed cells | 0.0 – 1.0                                |
| g_cond      | Condensate-driven ISG cytokine amplification     | Dimensionless coupling                   |
| R, K_R      | Counter-regulatory mediator & suppression const. | Dynamic state; K_R fitted                |
| L, η, v     | Latent pool, reactivation & seeding rates        | Dynamic state; η, v fitted               |

All new symbols are additive; setting φ_variant = 1, ε_drug = 0, π = 0, g_cond = 0, R ≡ 0, v = 0 recovers the core model exactly.

---

## 6. Theoretical Hypotheses & Falsification Criteria (Advanced)

**Hypothesis H5 (Cytokine Balancing)**  
The regulatory arm R suppresses peak pro-inflammatory cytokine C and alters safety-boundary trade-offs.  

*Falsified if:* profile rankings remain unchanged across sweeps of K_R and σ_R.

**Hypothesis H6 (Drug Selectivity)**  
Host-targeted entry inhibitors reduce the required TIP dose for extinction only when wild-type entry suppression exceeds TIP entry suppression.  

*Falsified if:* required TIP dose rises monotonically for all ε_drug > 0 under symmetric reduction.

**Falsifier F6 (Latent-Pool Accessibility)**  
If realistic latent parameters (v, η) render sterilising extinction unreachable, therapeutic outcomes must be reclassified from “viral clearance” to “long-term viral suppression.”

These criteria extend the core falsifiers F1–F4 without contradiction. Productive extinction (the quantity scored by the existing engine) remains the primary numerical target; sterilising extinction is an additional, stricter endpoint once the latent compartment is implemented.

---

## Consistency Notes

- Critical entry-blocking fractions quoted in the framework (88.8 % baseline, up to 97.2 % for 4-fold acceleration) are consistent in magnitude with the core analytic result 1 − 1/R₀ ≈ 0.885–0.888 (R₀ ≈ 8.7–8.9).  
- The Hill-function dose multiples for 89 % inhibition are pure arithmetic and do not conflict with any core parameter.  
- The primed-cell protection factor π generalises the core primed-cell postulate (P2); π = 0 recovers the core equations exactly.  
- Latent dynamics do not alter the productive-extinction criterion already implemented (V_W < V_ext).  

No overlapping or contradictory numerical claims have been retained. All statements above are structural or directly derivable from the Hill function / next-generation-matrix structure already verified in the core document.

