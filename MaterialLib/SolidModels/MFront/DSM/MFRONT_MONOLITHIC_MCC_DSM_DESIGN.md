# Monolithic MCC + DSM double-structure MFront behaviour — DESIGN

**Branch:** `dsm_mfront_maxwell` (off `dsm_mfront_hierarchical`, worktree
`~/git/ogs-worktrees/dsm_mfront_maxwell_wt`).
**Status:** DESIGN ONLY (no `.mfront` written yet). Agent B1, 2026-06-08.

This document specifies a **single** `@DSL ImplicitGenericBehaviour` that solves
Modified-Cam-Clay plasticity **and** the DSM micro/macro double-structure inside
**one** local Newton — every internal state variable in **one coupled implicit
residual vector**, populated and solved **once per Gauss point per global
iteration**, with the **integrable Maxwell web** (`mu_lR_mech`) as the
cross-coupling that makes the local tangent symmetric.

It does NOT make scientific decisions. The physics is taken verbatim from
(a) Beese's local Newton in `THM_DSM_Richards.nb` and (b) the already-derived,
GP-verified native web in
`dsm_native_pdisj_maxwell/.../ConstitutiveRelations/PotentialExchange.h`
(+ its `DSM/*.md` design docs). Every residual entry and Jacobian block below is
cross-referenced to **both** sources so the two can be checked against each other.
Anything genuinely ambiguous between the two is flagged **[ASSUMPTION]** and
proceeds with the documented default; the numerical/scientific call remains
Vinay's.

Tags: **[D]** derived/established in a cited source · **[ASSUMPTION]** documented
default where the two sources do not pin it · **[PRED]** predicted, not verified.

---

## 0. The two source skeletons, side by side

### 0.1 Beese local Newton (`THM_DSM_Richards.nb`, "Build Residual" cell) [D]
The notebook assembles, per GP, the residual vector

```
R_g = Flatten[{ R_epsth, R_phiM, R_nl, R_rhoSR, R_rholR, R_epssw, R_epsp, R_pc }]
```

solved for the unknown set (one residual entry ⇄ one unknown):

| # | unknown (Beese) | residual `R_*` (Beese, decoded from the nb) | meaning |
|---|---|---|---|
| 1 | `epsth` (θ-vol strain) | `epsth − (epsthn + Δt·c[Liquid][TempRate]·𝒫[Solid][α_th])` | thermal volumetric strain rate eq. |
| 2 | `phiM` (macro porosity) | `phiM − (phiMn + Δt·c[Solid][History][phiMdot])` | macro-porosity evolution (solid mass bal.) |
| 3 | `nl` (micro water content) | `storagel + phiPenalty·⟨−(phiM−phiMmin)⟩²` ; `storagel := ρ̇_l − α_M·ρ̂_l − ρ_l·ε̇` | **micro storage / exchange** + macro-porosity floor penalty |
| 4 | `rhoSR` (solid real density) | `rhoSR − (rhoSRn + Δt·c[Solid][History][rhoSRdot])` | solid-grain EOS rate |
| 5 | `rholR` (liquid real density) | `rholR − c[Liquid][History][rholR]` | micro-liquid EOS (algebraic) |
| 6 | `epssw` (swelling strain) | `epssw − (epsswn + Δepssw)` | swelling-strain increment |
| 7 | `epsp` (plastic strain, Voigt) | `Voigt[epsp − (epspn + Δepsp)]` | MCC plastic-flow increment |
| 8 | `sigpc` (preconsolidation) | `sigpc − (sigpcn + Δpc)` | MCC hardening increment |

Scaling: `Φ := Sqrt[R_g·R_g] + c[Solid][History][Φ]`; objective `𝕈g = R_g` (the
`resScaleVec` per-row scaling is present but commented out in the shipped nb).
This is the **template** the MFront `@Integrator` reproduces — same eight physical
unknowns, re-expressed in MFront/TFEL idiom.

### 0.2 Native web (`PotentialExchange.h`, `dsm_native_pdisj_maxwell`) [D]
The GP-level physics already implemented and verified:

- `computeVanDerWaalsMicroPotential` → `mu_lR_vdw = −Π/ρ_lR` [J/kg],
  `Π = (A·Sa³/6π)·nS³ρ_SR³/(n_l³ρ_lR)` plus the lumped augmentation
  `K·exp(−h/λ)` (ADDITIVE, `+=`, never overwrite — §4.1 CLAUDE.md), with analytic
  `dmu_lR_dnl`, `dmu_lR_drho_lR`, `d2mu_lR_dnl2`.
- `computeIntegrableMechanicalMicroPotential` (the **web**, spec item 2):
  ```
  mu_lR_mech = −[ (Π + n_l·Π')·eps_v + 0.5·b·K_drained·eps_v² ] / ρ_lR   [J/kg]
  mu_lR      = mu_lR_vdw + mu_lR_mech     (ADDITIVE)
  ```
  with analytic `dmu_lR_mech/{deps_v, dn_l, drho_lR}` (L485–513). Maxwell identity
  it satisfies (L463–471): `d sigma_sw,m/d n_l = n_S·ρ_lR·d mu_lR/d eps_v` on the
  drained line `p_conf = −K_drained·eps_v`.
- `computePotentialDrivenMassExchange` → `ρ̂_l = α_M·(mu_LR − mu_lR)`, `ρ̂_L = −ρ̂_l`.
- Eigenstress (REV mean): `sigma_sw,m = −φ_m·p_film = −n_S·n_l·(Π − b·p_conf)`,
  `φ_m = n_S·n_l`, `n_S = 1 − φ_M`, `p_conf = −tr(σ')/3`.
- REV referencing corrigendum (`MAXWELL_CONJUGATE_REV_REFERENCING.md`): the
  *strain-view* conjugate divides by `ρ_lR·n_S`; the **integrable** partner
  `mu_lR_mech` (the form this design uses) carries **no** explicit `(1−φ_M)` —
  the contact-area `(1−φ_M)` cancels the per-REV-mass referencing, so it divides
  by the intrinsic `ρ_lR` (equipresence note, `PotentialExchange.h` L351–358).

### 0.3 Term-by-term correspondence (the check Vinay asked for)

| physical content | Beese nb symbol | native web symbol | this MFront `@StateVariable` |
|---|---|---|---|
| micro water content | `nl` | `n_l` | `n_l` |
| micro liquid real density | `rholR` | `rho_lR` | `rho_lR` |
| solid real density | `rhoSR` | `rho_SR` (param here, const grain) | `rho_SR` (aux; see §4) |
| macro porosity | `phiM` | `phi_M` (= φ − φ_m hierarchical) | `phi_M` |
| swelling strain | `epssw` | `epsilon_sw` (= `swelling_slope·Δφ_m`) | `epsilon_sw` |
| plastic strain | `epsp` | (MCC, native carries via FEM) | `eps_el` (elastic; plastic implied) |
| preconsolidation | `sigpc` | `pc` | `rpc` (`pc = rpc·E`) |
| micro chem. potential | `mu_lR` | `mu_lR = mu_lR_vdw + mu_lR_mech` | (computed, aux `mu_lR_value`) |
| disjoining pressure | `p_L_m` | `Π = −ρ_lR·mu_lR_vdw` | (computed) |
| **mech. Maxwell partner** | (absent — Beese has no `mu_lR_mech`) | `mu_lR_mech` | (computed; §5) |
| mass exchange | `ρ̂_l` (`storagel`) | `rho_l_hat = α_M(μ_LR−μ_lR)` | (computed, aux `rho_l_hat_value`) |
| thermal vol strain | `epsth` | (isothermal native; absent) | **omitted** — see §3 note |

**Two structural deltas between the sources, both documented:**
1. Beese carries `epsth` (thermal) and a live `rhoSR` EOS; the native web is
   isothermal with constant grain. This MFront mirrors the **native** scope
   (isothermal, const grain) — `epsth` dropped, `rho_SR` a parameter — because
   that is the verified web and the ANCHORS benchmark switches both channels off
   (`MAXWELL_CONJUGATE_IMPLEMENTATION.md` §6.1). **[ASSUMPTION]** isothermal +
   rigid grain; revisit (B1.5) only for a thermal / `α_B<1` benchmark.
2. Beese has **no** `mu_lR_mech` (his `μ_lR` is adsorption-only). The native web
   ADDS the integrable partner. This MFront includes it — that **is** the point
   of the `maxwell` branch. So this behaviour = Beese's monolithic Newton
   **structure** + the native web's `mu_lR_mech` **physics**.

---

## 1. DSL, gradients, thermodynamic forces

```
@DSL ImplicitGenericBehaviour;
@Behaviour RichardsMechanicsDSMMicroMacroMaxwell_MCC;
@Algorithm NewtonRaphson;     // (or LevenbergMarquardt fallback — see §8)
@Theta 1;  @Epsilon 1e-14;  @MaximumNumberOfIterations 200;

@Gradient StrainStensor eps_to;   eps_to.setGlossaryName("Strain");
@ThermodynamicForce StressStensor sig;  sig.setGlossaryName("Stress");
@Gradient real p_LR;  p_LR.setEntryName("LiquidPressure");
@ThermodynamicForce real S_L;  S_L.setEntryName("Saturation");
```

Same external coupling as the existing `_MCC` bridge (strain+pressure in,
stress+saturation out), so the RM process driver is unchanged. The DIFFERENCE
is internal: the micro state is promoted from an `@InitializeLocalVariables`
hand-rolled inner solver (current bridge) to **first-class `@StateVariable`s in
the one MFront implicit system** (this design).

---

## 2. The monolithic `@StateVariable` list (the one Newton's unknowns)

The local Newton solves for the increments of all of these **together** (TFEL
appends `Δ` and the implicit system is `f<var> = 0`):

```
// ── mechanical (MCC) ──────────────────────────────────────────────
@StateVariable StrainStensor eps_el;   eps_el.setGlossaryName("ElasticStrain");   // (6 / 4 comp)
@StateVariable real          Lam_p;    Lam_p.setGlossaryName("EquivalentPlasticStrain");
@StateVariable strain        rpc;      // scaled preconsolidation, pc = rpc*young

// ── DSM micro/macro double structure ──────────────────────────────
@StateVariable real n_l;       n_l.setEntryName("n_l");        // micro water content
@StateVariable real rho_lR;    rho_lR.setEntryName("rho_lR");  // micro liquid real density
@StateVariable real phi_M;     phi_M.setEntryName("phi_M");    // macro porosity
@StateVariable real epsilon_sw; epsilon_sw.setEntryName("epsilon_sw"); // swelling strain (scalar, isotropic)
```

Unknown-count: `eps_el` (Stensor, 4 or 6) + `Lam_p` (1) + `rpc` (1) + `n_l` (1)
+ `rho_lR` (1) + `phi_M` (1) + `epsilon_sw` (1) = **N_mech + 5 scalars**.

Mapping to Beese's eight (§0.1):
- `eps_el` ⇄ Beese `epsp` (complementary: `eps_p = eps_to − eps_el`; MFront's
  StandardElasticity idiom carries the *elastic* strain, Beese carries plastic —
  same information).
- `rpc` ⇄ `sigpc`. `Lam_p` ⇄ the MCC plastic multiplier `Λ` (Beese folds flow
  into `Δepsp`; MFront's MCC idiom carries `Λ_p` + flow direction — see §6).
- `n_l, rho_lR, phi_M, epsilon_sw` ⇄ Beese `nl, rholR, phiM, epssw`.
- Beese `epsth, rhoSR` → **dropped** (isothermal, const grain; §0.3 note 1).

**Auxiliaries (reported, not solved):** `phi_m`, `phi`, `n_S`, `n_L`, `rho_LR`,
`omega_l`, `mu_lR`, `rho_l_hat`, `delta_epsilon_sw`, `sigma_S`, `pc`,
`PlasticVolumetricStrain`, `VolumeRatio` — exactly the `@AuxiliaryStateVariable`
set the current `_MCC` bridge already exposes (L43–92), so post-processing and
the parity runner see the same outputs.

---

## 3. The coupled residual vector (one entry per unknown)

Let `Δt`, `θ` per TFEL. `tr(Δeps_to) = Δeps_v` (volumetric strain increment).
`p = p_LR + Δp_LR`. All residuals are written so that the **converged** state
reproduces, term for term, the native web AND Beese's `R_*`.

### r1 — `feps_el` (mechanical, MCC return map) [D: Beese R_epsp + native FEM MCC]
```
feps_el = Δeps_el + Δeps_pl − Δeps_to        (elastic ⇄ Beese R_epsp)
```
with `Δeps_pl = ΔLam_p · n_flow`, `n_flow` the MCC normalized flow direction
(identical algebra to the current bridge `@Integrator` L1385–1390). The swelling
eigenstress enters the **stress** (§7), not this kinematic residual — swelling is
a stress, not a strain (`MAXWELL_PAIR_RESTORATION.md` §1: β_sw eigenstrain
retired). `epsilon_sw` is carried only to report the volumetric swelling and to
drive the eigenstress; see r7.

### r2 — `fLam_p` (MCC yield) [D: native/current bridge L1391]
```
fLam_p = f / fchar ,   f = q² + M²·p_mcc·(p_mcc − pc) ,   p_mcc = −tr(σ)/3 + pamb
fchar = pc_char · young
```
KKT handled by `@AdditionalConvergenceChecks` (elastic switch if `ΔLam_p < 0`),
as in the current bridge (L1454–1467).

### r3 — `frpc` (MCC hardening) [D: Beese R_pc; native L1392]
```
frpc = Δrpc + Δeps_pl_v · θ_h · (rpc_new − rpc_min) ,  θ_h = vr/(la−ka)
```
(`vr` = volume ratio aux; `Δeps_pl_v = tr(Δeps_pl)`.)

### r4 — `fn_l` (micro storage / mass exchange) [D: Beese R_nl `storagel`; native exchange]
This is the heart of the coupling. Beese: `storagel := ρ̇_l − α_M·ρ̂_l − ρ_l·ε̇`.
Discrete, REV mass form (matching the current bridge's `residual_1`, L712):
```
fn_l = φ_m·ρ_lR − (φ_m·ρ_lR)|_n − Δt·ρ̂_l − φ_m·ρ_lR·Δeps_v
```
with
```
φ_m   = n_S·n_l = (1−φ_M)·n_l            (hierarchical; native φ_m = n_S·n_l)
ρ̂_l   = α_M·(μ_LR − μ_lR)                (computePotentialDrivenMassExchange)
μ_LR  = 0 (p≥−tol) | p/ρ_LR  (else)      (computeYoungLaplaceMacroPotential)
μ_lR  = μ_lR_vdw(n_l,ρ_lR,n_S) + μ_lR_mech   (§5; ADDITIVE — never overwrite)
```
**[ASSUMPTION]** Beese's macro-porosity floor penalty
`phiPenalty·⟨−(phiM−phiMmin)⟩²` is reproduced as a **clamp/penalty on r6**
(`phi_M ≥ phi_M_min`), not folded into r4, to keep r4 a clean mass balance. The
current bridge enforces the same floor by clamping `phi_M = clamp(φ−φ_m, 0, φ)`
(L1347). Default: keep the clamp; expose `phiPenalty` only if a benchmark needs
the smooth penalty. Vinay's call if the penalty must be smooth.

### r5 — `frho_lR` (micro-liquid EOS) [D: Beese R_rholR; native `rho_micro_eos`]
```
frho_lR = rho_lR − rho_micro_eos(n_l, rho_lR)
rho_micro_eos = rho_l0·exp(−a·ω^b) + rho_LR_ref ,  ω = n_l·ρ_lR/(n_S·ρ_SR)
```
(identical to current bridge L390–396; `a`,`b` are `MicroLiquidDensityA/B`.)

### r6 — `fphi_M` (macro-porosity evolution) [D: Beese R_phiM]
Beese: `phiM − (phiMn + Δt·phiMdot)`. In the native **hierarchical** split,
`phi_M` is not independent — it is `φ_total − φ_m`, with `φ_total` advected by the
volumetric strain (current bridge L324–338, `phi_total_trial`). Two equivalent
encodings; this design uses the algebraic hierarchical closure as the residual:
```
fphi_M = phi_M − (phi_total_trial − φ_m) ,   φ_m = (1−φ_M)·n_l   (implicit in φ_M)
phi_total_trial = clamp( (phi_total_prev + Δeps_v)/(1+Δeps_v), n_l_min, 1−ε )
```
**[ASSUMPTION]** algebraic hierarchical closure (matches the verified native +
current-bridge porosity split) rather than Beese's rate form `phiMdot`. They
coincide when `phiMdot` is the solid-mass-balance rate at `α_B=1`, const grain
(`MAXWELL_CONJUGATE_IMPLEMENTATION.md` §6.1). The floor `phi_M ≥ phi_M_min`
(Beese `phiMmin`) is applied here.

### r7 — `fepsilon_sw` (swelling strain) [D: Beese R_epssw; current bridge L1349]
```
fepsilon_sw = epsilon_sw − (epsilon_sw|_n + Δepsilon_sw) ,
Δepsilon_sw = swelling_slope · (φ_m − φ_m|_n)
```
(scalar isotropic; drives the spherical eigenstress in §7.)

**Residual ordering** (Jacobian block layout, §4):
`{ feps_el, fLam_p, frpc, fn_l, frho_lR, fphi_M, fepsilon_sw }` — mechanical block
first (as in the current bridge / MCC idiom), DSM block second. This is Beese's
`R_g` order with `epsth/rhoSR` removed and the MCC trio (`feps_el, fLam_p, frpc`)
expanded from his single `epsp/pc` pair to the MFront return-map idiom.

---

## 4. Jacobian block structure (the symmetric web)

The local tangent `J = ∂f/∂Δ(unknowns)`. Write the unknown vector as
`X = [eps_el | Lam_p | rpc || n_l | rho_lR | phi_M | epsilon_sw]`.
Block form (M = mechanical 3×3 super-block, D = DSM 4×4 super-block, C = couplings):

```
        eps_el  Lam_p  rpc  | n_l    rho_lR  phi_M  eps_sw
feps_el [ Mee    Mep   Mer  | 0       0       0      0    ]
fLam_p  [ Mpe    0     Mpr  | 0       0       0      0    ]
frpc    [ Mre    Mrp   Mrr  | 0       0       0      0    ]
        --------------------+----------------------------
fn_l    [ C_ne   0     0    | D_nn    D_nr    D_np   0    ]   ← C_ne = ∂fn_l/∂eps  (THE WEB)
frho_lR [ 0      0     0    | D_rn    D_rr    0      0    ]
fphi_M  [ C_Me   0     0    | D_Mn    0       D_MM   0    ]
feps_sw [ 0      0     0    | C_sn    0       C_sM   1    ]
```

- **M-block (mechanical):** identical to the current bridge `@Integrator` analytic
  Jacobian (L1394–1414): `∂feps_el/∂Δeps_el`, `∂fLam_p/∂Δeps_el`, the MCC `∂n/∂…`.
- **D-block (DSM):** `D_nn, D_nr, D_rn, D_rr` are exactly the current bridge inner
  Newton `j11,j12,j21,j22` (L1047–1052), now promoted to the global MFront
  Jacobian instead of an inner hand-solve. `D_np = ∂fn_l/∂φ_M` and
  `D_Mn = ∂fphi_M/∂n_l` come from `φ_m = (1−φ_M)·n_l`.
- **C_ne = ∂fn_l/∂Δeps (THE INTEGRABLE WEB COUPLING).** Through `μ_lR_mech`:
  ```
  ∂fn_l/∂eps_v = −Δt·∂ρ̂_l/∂eps_v = −Δt·α_M·(−∂μ_lR/∂eps_v)
               = +Δt·α_M·∂μ_lR_mech/∂eps_v
  ∂μ_lR_mech/∂eps_v = −[ (Π + n_l·Π') + b·K_drained·eps_v ] / ρ_lR   (native L505)
  ```
  PLUS the explicit `−φ_m·ρ_lR·Δeps_v` storage term's `∂/∂eps_v = −φ_m·ρ_lR`
  (current bridge `∂r1/∂eps_v`, L1197). As a Kelvin/Voigt row: `C_ne ∝ mᵀ`
  (`m` = identity Voigt), since only the volumetric part `eps_v = mᵀeps` enters.
- **Mechanical ← DSM coupling (the OTHER half of the pair).** The swelling
  eigenstress `sigma_sw,m = −φ_m·(Π − b·p_conf)` enters the **stress** (§7), so it
  shows up in `feps_el`/`fLam_p` through `σ`. Its `n_l`-derivative
  `∂sigma_sw,m/∂n_l = S₁ = −n_S·(Π + n_l·Π')` is the **transpose partner** of
  `C_ne` (native L240–317, REV-referencing doc §3). **The integrability the
  branch is named for:** `C_ne` (lower-left, `∂μ_lR/∂eps`) and the eigenstress
  upper-right block (`∂σ/∂n_l`) are a transpose pair (both ∝ `S₁`) ⇒ the local
  tangent is **symmetric** ⇒ derives from one `Ψ`. With the integrable
  `mu_lR_mech` (quadratic-in-eps_v form), this holds **on and off** the drained
  line, not just at `p_conf = −K_drained·eps_v` (native L463–471).

**Monolithic payoff over the current bridge:** the current `_MCC` bridge solves
the DSM micro state in a **hand-rolled inner Newton** inside
`@InitializeLocalVariables` (L427–1171, ~750 lines: bracketing, line-search,
explicit fallbacks) and then injects only `dn_l/deps_v` and `dn_l/dp` into the
tangent (L1497–1519). This design **deletes that inner solver**: `n_l, rho_lR,
phi_M, epsilon_sw` become MFront implicit unknowns, TFEL's outer Newton solves
the whole `(mechanical, micro, macro)` system at once, and the Maxwell coupling
`C_ne` is a first-class Jacobian block — populated once, solved once
(spec: "the ENTIRE local Newton populated and solved ONCE"). This mirrors Beese's
single `FindRoot` over `R_g` exactly.

---

## 5. The Maxwell web partner `mu_lR_mech` (verbatim from the native web)

Computed each residual evaluation (NOT a stored state):
```
mu_lR     = mu_lR_vdw + mu_lR_mech                                   // += , J/kg
Π         = −ρ_lR·mu_lR_vdw                                          // disjoining, Pa
Π'        = dΠ/dn_l   (from native dmu_lR_dnl: Π' = −ρ_lR·dmu_lR_vdw_dnl)
Π''       = d²Π/dn_l² (from native d2mu_lR_dnl2)
p_conf    = −tr(σ')/3                                                // confining, Pa
K_drained = young/(3(1−2ν))                                          // already in bridge as K
b         = biot_coefficient                                          // [ASSUMPTION] param; ANCHORS b=1
mu_lR_mech = −[ (Π + n_l·Π')·eps_v + 0.5·b·K_drained·eps_v² ] / ρ_lR  // native L503
```
Derivatives `∂mu_lR_mech/{∂eps_v, ∂n_l, ∂rho_lR}` taken verbatim from native
L505–511. **`+=` discipline (CLAUDE.md §4.1):** `mu_lR` accumulates
`mu_lR_vdw` then `mu_lR_mech`; replacement `=` is forbidden (overwriting-vdW-base
incident). **Units annotated on every line (§4.2).**

Sign/gate decisions are INHERITED, not re-made (all Vinay's, 2026-06-02,
`MAXWELL_CONJUGATE_IMPLEMENTATION.md` §6):
- load EXPELS micro water; OGS effective stress tension-positive;
- the integrable form (this design) uses the **smooth quadratic** energy, so the
  sharp `⟨p_conf−Π⟩₊` gate of the strain-view helper is **superseded** — the web
  is `C¹` by construction (native L447–476). **[ASSUMPTION]** use the integrable
  partner (no sharp gate), per the `maxwell` branch's stated spec; if Vinay wants
  the sharp-gate film-pressure variant instead, that is `computeFilmPressure…`
  (native L369) — a different branch decision.
- `b = biot_coefficient`: ANCHORS sets `b=1` (rigid grain, §6.1). Carried as a
  `@MaterialProperty`, cited to the PRJ (§12 provenance), NOT hard-coded.

---

## 6. MCC yield + hardening (unchanged from the verified bridge)

Yield `f = q² + M²·p(p−pc)`, associated flow, `pc` hardening
`ṗc = pc·(ε̇_p_v)/(λ−κ)·v_r`. Algebra, normalization, and the elastic/plastic
switch are taken **verbatim** from `ModCamClay_semiExpl_constE` /
`RichardsMechanicsDSMMicroMacroBridge_MCC` `@Integrator` (L1358–1467). No new MCC
physics — the only change is that `pc`'s residual `frpc` now lives in the **same**
Newton as the DSM block, so plastic compaction (`Δeps_p_v`) couples to micro
expulsion through `eps_v` in `C_ne` within one solve (Beese's `R_pc` and `R_nl`
are siblings in his `R_g`).

---

## 7. Stress update — eigenstress carries the swelling (and the Maxwell LEFT half)

```
@ComputeFinalThermodynamicForces (schematic):
σ = σ₀ + ∂σ/∂eps_el · Δeps_el  +  sigma_sw,m · I
sigma_sw,m = −φ_m·(Π − b·p_conf) = −φ_m·p_film          // native eigenstress, Pa
```
The current bridge writes `σ = … − K·Δepsilon_sw·I` (L1419), i.e. swelling via the
`epsilon_sw` increment × bulk modulus — the **eigenstrain-flavoured** encoding.
This design keeps that encoding for `epsilon_sw` consistency BUT the **Maxwell
left half** is `∂σ_sw,m/∂n_l = S₁`; whichever encoding is used, the stress's
`n_l`-derivative must equal `S₁` so it transposes `C_ne` (§4). **[ASSUMPTION]**
keep the `−K·Δepsilon_sw·I` form (verified, parity-tested) and verify
`∂σ/∂n_l == S₁` in the GP Maxwell-symmetry test rather than rewriting the stress
update; if they differ, that is a finding for Vinay, not a silent re-encode.

---

## 8. Integration scheme

- `@DSL ImplicitGenericBehaviour`, `@Algorithm NewtonRaphson`, analytic Jacobian
  (all blocks §4 available in closed form from the two sources). **[PRED]** the
  full analytic Jacobian gives quadratic local convergence (not benchmarked).
- Fallback: `@Algorithm LevenbergMarquardt` or a numerical-Jacobian build
  (`@CompareToNumericalJacobian` during bring-up) to validate the analytic blocks
  — the standard MFront verification path; the current bridge's bracketing /
  line-search robustness (L468–584) is no longer needed once the micro state is in
  the outer Newton, but `@Predictor` / bounds on `n_l∈(0,φ]`, `rho_lR>ρ_LR_ref`,
  `phi_M∈[phi_M_min,φ]` must be set (`@Bounds`/`@PhysicalBounds`).
- `@Theta 1` (fully implicit, matches Beese's backward-Euler `… n + Δt·(…)`).
- Tolerances derived from problem scale (CLAUDE.md §1.2/§3) — NOT raw literals;
  the pressure-scale floor incident (≥1 Pa at 100 MPa) applies to the
  `p_LR`-coupled rows. **[ASSUMPTION]** reuse the current bridge's `@Epsilon 1e-14`
  on the scaled residual; confirm against the 100 MPa suction scale before any run.

---

## 9. What this design deliberately does NOT decide (Vinay's calls)

1. Material/parameter VALUES — none here; all live in PRJ/MTEST with §12
   provenance + §1.1 citations. Structure is parameter-free.
2. Whether to keep the `−K·Δepsilon_sw·I` eigenstrain encoding or switch to a
   direct `−φ_m·p_film·I` eigenstress (§7) — verify-then-ask.
3. Sharp-gate film-pressure variant vs integrable web (§5) — the branch name says
   integrable; flagged in case.
4. Re-introducing `epsth` (thermal) / live `rho_SR` (Beese's full set) for a
   thermal or `α_B<1` benchmark (B1.5, §0.3) — out of scope until such a case.
5. The macro-porosity floor as hard clamp vs smooth penalty (§3 r6).
6. All test EXPECTED values and tolerances (§3 CLAUDE.md) — TODO(user).

## 10. Verification anchors (structure only; no expected values)

Per CLAUDE.md §3 — proposed test STRUCTURE, expected values TODO(user):
- **Maxwell-symmetry GP test** [anchor: derived identity]: FD-confirm
  `n_S·ρ_lR·∂μ_lR/∂eps_v == ∂σ_sw,m/∂n_l` at a sample state (native L463; REV doc
  §3). Calibration-free.
- **Below-gate / free-swelling regression** [anchor: approved baseline, cite the
  MS LE ANCHORS reference VTU]: reproduce the current `_MCC` bridge result where
  the web term is inert.
- **Parity vs native** [anchor: approved baseline]: add a pair to
  `scripts/run_dsm_parity.py` PARITY_SUITES; assert MAE vs
  `dsm_native_pdisj_maxwell`.
- **Numerical-Jacobian check** [anchor: integrability]: `@CompareToNumericalJacobian`
  on every §4 block during bring-up.

## 10a. Why monolithic closes the audit's central violation (and the ONE open call it does NOT)

`EQUIPRESENCE_AUDIT_2026-06-06.md` finds the native web's **central** defect is
that `mu_lR_mech` is added only at the **2 global** exchange-assembly sites, not at
the **3 micro local-solve** sites — so each GP converges `n_l` against a
gate-closed `mu_lR` while the macro residual it feeds uses the gate-open `mu_lR`:
**two values of one potential per GP** ("two functions sharing a name", audit
§3.1). A correct `mu_lR(p_film)` must appear at *every* balance, including the one
that determines its own argument `n_l`.

**A monolithic MFront Newton makes this violation structurally impossible.** There
is no separate inner micro solve to fall out of sync: `n_l, rho_lR, phi_M` are
unknowns of the **same** implicit system whose residuals `fn_l, frho_lR` and whose
mechanical residuals all read **one** `mu_lR = mu_lR_vdw + mu_lR_mech`, evaluated
once per residual call (§5). This is precisely the audit's recommended fix
("one evaluator … the *same* `mu_lR(p_film)` returned to the micro residual, the
eigenstress driver, and the macro exchange", §3.3) realised by construction rather
than by re-routing three native call sites. It also retires the audit's
tangent-consistency rows that arise from the native local-vs-global Jacobian
mismatch (the `computeImplicitNlDpL` linearizing the wrong residual; the dropped
`p_L` channel), because the local tangent IS the global Jacobian block (§4).

**The ONE thing monolithic does NOT decide — and must stay Vinay's call: the
sign.** The audit verifies two independent ways (rate/residual + steady-state
`rho_l_hat=0`) that with the **settled, correct** swelling sign
`sigma_sw,m = −n_S·n_l·Π < 0`, integrability **forces** `S₁ > 0`, hence under
compression (`eps_v < 0`) `mu_lR_mech < 0` and the micro phase **gains** water —
the **opposite** of the documented "load EXPELS micro water" intent
(`MAXWELL_PAIR_RESTORATION.md §3`). This is not a wiring bug a monolithic solve
fixes; it is a statement that *adding `sigma_eff` with the integrability-consistent
sign produces inflow under load*. The two cannot both hold: either (a) "load
expels" → needs `S₁ < 0`, an anti-swelling term that **breaks the verified
swelling closure**, or (b) integrability-consistent inflow → the **documentation's
expulsion language** (and the gate direction) is what to revise. **[ASSUMPTION]**
This design wires the integrability-consistent sign (option b) — it preserves the
verified eigenstress and the one-`Ψ` symmetry the branch is named for — and flags
the contradiction here. The physics call (a vs b) is **Vinay's**; the GP
Maxwell-symmetry test (§10) plus a stress-loading expulsion test lock whichever he
chooses. Until he rules, no benchmark should be trusted on a gate-open / high-load
path (audit §2 sign-convention row).

## 11. Cross-references
- Beese: `~/git/GitHub/matModels_Beese/src/TPM/THM_DSM_Richards.nb` ("Build
  Residual" / "Scale Residual" cells; AceGen export `TRM_camclay_gp.c`,
  `hg[33]` history).
- Native web: `dsm_native_pdisj_maxwell/.../ConstitutiveRelations/PotentialExchange.h`
  (`computeIntegrableMechanicalMicroPotential` L485; vdW L106; exchange L528) and
  `DSM/{MAXWELL_CONJUGATE_IMPLEMENTATION,MAXWELL_PAIR_RESTORATION,
  MAXWELL_CONJUGATE_REV_REFERENCING,EQUIPRESENCE_AUDIT_2026-06-06}.md`.
- Current MFront bridge being superseded: `RichardsMechanicsDSMMicroMacroBridge_MCC.mfront`
  (the hand-rolled inner micro solver this design folds into the outer Newton).
- Guardrails: repo `CLAUDE.md` §1 (literals), §4 (`+=`/units), §12 (PRJ provenance).
