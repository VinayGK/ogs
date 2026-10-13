# Audit: appendix "The level form: what is derived, what is posited, what is missing"

Adversarial, read-only audit, 2026-10-13. Object: the LaTeX appendix `sec:BGR_app_derivation`.
Code: tag `submission-config-anchors-T4.1-T4.2-2026-10-09` = commit `8bbcace35c`. All line numbers are
read from `git show 8bbcace35c:<path>`. Symbolic checks: `audit_level_form_sympy.py` (sympy 1.14),
output in `sympy_output.txt` (17/17 OK). Abbreviations: FEM = `ProcessLib/RichardsMechanics/RichardsMechanicsFEM-impl.h`,
PE = `.../ConstitutiveRelations/PotentialExchange.h`, PEP = `.../PotentialExchangeParameters.h`,
CRP = `.../CreateRichardsMechanicsProcess.cpp`.

**Line numbers.** No cited range is shifted: at the tag, every cited range contains the cited code.
The minor citation problems are listed below at 1.5, 1.12 and 1.16; none is an offset.

## (1) Statements labelled *implemented*

1.1 **HOLDS**. The level is at FEM l.3571–3572: `double L_curr = -n_S * n_l * (Pi_curr_film - b_film * p_conf_film);`, so L = −n_S n_l [Π + b σ'_mean] with p_conf = −σ'_mean. Note: this line is the *lagged* form, with `p_conf_film` taken from the previous Newton evaluation (comment l.3553–3555). The closed form at l.3573–3586 overwrites it whenever m̂ is finite.

1.2 **HOLDS**. The closed form is at FEM l.3573–3586. Comment l.3581: `L (1 + c) = -n_S n_l Pi - c (m_hat + s_prev - L_prev)`. Code l.3582–3585: `c_lv = n_S * n_l * b_film; L_curr = (-n_S*n_l*Pi_curr_film - c_lv*(level_m_hat + sigma_sw_prev_mean - L_prev)) / (1.0 + c_lv);`. Sympy: solving L = −φ_m[Π + b(m̂ + s_prev + L − L_prev)] for L gives this expression exactly.

1.3 **HOLDS**. FEM l.3591: `delta_sigma_sw.noalias() += (L_curr - L_prev) * identity2_film;`.

1.4 **HOLDS**. For the previous level, FEM l.3566–3568 is `L_prev = std::isfinite(level_prev_used) ? level_prev_used`, and the fallback at l.3569–3570 is `: -n_S_prev_level * n_l_prev * (Pi_prev_film - b_film * p_conf_prev_film)`. The stored level is read from `PrevState<SwellingLevelUsed>` (l.3983–3987) and written at l.4050–4053. Note: a stored value of exactly 0.0 counts as "not set" (l.3985).

1.5 **HOLDS**. L⁽⁰⁾ comes from the initial-state evaluation, FEM l.9126–9143: `initial_level_reference = ... **PrevState<SwellingLevelUsed> == 0.0`, followed by `**PrevState<SwellingLevelUsed> = *SwellingLevelUsed(current)` (l.9140–9143). Citation gap: the claim that it "books no stress" rests on l.4039–4049 (`if (level_form_state && initial_reference) { sigma_sw = *sigma_sw_prev; }`), and that range is not cited. The telescoping σ_sw⁽ⁿ⁾ = L⁽ⁿ⁾ − L⁽⁰⁾ also needs σ_sw = 0 at t0. The appendix does not state this.

1.6 **HOLDS**. A rejected step reuses the level of the last accepted step, because the previous level is read only from `PrevState<SwellingLevelUsed>` (l.3985–3986), and the current evaluation writes only to `state_current` (l.4052).

1.7 **HOLDS**. The current level takes K at the current total porosity: FEM l.3250–3251 `K_aug_sw = effectiveAugmentationPrefactor(potential_exchange_params, total_porosity)`, called with `PorosityData.phi` (l.4033–4035). The comment at l.418 gives φ = φ_M + (1−φ_M) n_l.

1.8 **HOLDS**. The fallback at FEM l.3258–3262 evaluates `K_aug_sw_prev` at `total_porosity_prev` when `per_level_K_sw` is set. The level form forces this flag on (PEP l.781–785: `return p.swelling_stress_K_level || p.swelling_stress_form == SwellingStressForm::Level;`), and the previous porosity is supplied at FEM l.3949–3957.

1.9 **HOLDS**. PEP l.971–984 is `effectiveAugmentationPrefactor`. At l.978–981 it calls `getValueLogLinear(micro_solid_density_reference * (1.0 - phi))`. PEP l.64–86 then computes `K_l * std::exp(t * std::log(K_r / K_l))`, so ln K is linear in ρ_d between the knots and K is held flat outside the table.

1.10 **HOLDS**. μ_mech is at PE l.584–587: `A = Pi + n_l*dPi_dnl; Bmod = biot_b*K_drained; out.mu_lR_mech = -(A*eps_v + 0.5*Bmod*eps_v*eps_v)/rho_lR;`. Sympy confirms this equals Eq. BGR_dv_partner.

1.11 **HOLDS** (with qualification). The fold point is FEM l.1002: `out.mu_lR += mech.mu_lR_mech; // J/kg, additive`. It is "ungated" in the sense that no saturation or Π gate applies. It is reached only when `film_pressure_coupling` is set and p_conf is finite (l.811–815), and only when `film_strain_coupling == Off` (otherwise there is an early return at l.914–972). The level form requires both conditions at create time (CRP l.1107–1129), so for the level form the qualification is met.

1.12 **HOLDS** (minor citation). The m̂ prediction is m̂ = m_lag − s_lag + K(ε_v − ε_v,lag), coded at FEM l.4016–4018: `level_m_hat_sw = -p_conf_swelling - s_lag + K_d_sw * (variables.volumetric_strain - (v_lag - 1.0));`. s_lag is `SwellingLagStress` (l.4012–4013). The cited l.3992 is a comment line (`m_lag = -p_conf ...`) inside the comment block l.3990–3996.

1.13 **HOLDS**. The first evaluation uses the lagged level, FEM l.4003: `if (v_lag != 0.0 && std::isfinite(p_conf_swelling))`. With v_lag = 0 (unset), m̂ stays NaN and the code takes the lagged form.

1.14 **HOLDS**. "At convergence this is Eq. BGR_dv_closed": substituting m̂ = σ0 + K ε_v and s_prev − L_prev = −L⁽⁰⁾ into the code's form reproduces the equation exactly (sympy). Note: for a yielding skeleton an additional correction loop runs (FEM l.3784–3860, "NOT FOR PRODUCTION"). It does not run for a linear-elastic skeleton (l.3797–3798), and the appendix restricts itself to the linear case.

1.15 **HOLDS**. The exchange law is at PE l.985: `out.rho_l_hat = alpha_M * (mu_LR - mu_lR);`.

1.16 **HOLDS**. The KKT branch: FEM l.1949–1954 returns Interior unless `at_wall.f < 0`. At l.1987–1998 it sets `n_l = n_max`, `rejected_exchange = max(0, rho_l_hat - S_s)`, `exchange_received = S_s` and `booked_storage_rate = S_s`.

1.17 **CANNOT CHECK**. The submitted configuration sets `micro_ceiling_treatment = kkt`. The option exists (CRP l.93, l.951–953), but no submitted deck is in the repository.

1.18 **HOLDS**. The label: FEM l.7677–7678 reads `DIAGNOSTIC fix (b), level form: drain-feedback tangent (2026-10-04, JACOBIAN-ONLY; NOT FOR PRODUCTION)`. Only this block carries "NOT FOR PRODUCTION". The other blocks the level form adds or rescales (l.7635–7652, l.7735ff.) are labelled "DIAGNOSTIC fix (b)" or "JACOBIAN-ONLY" only.

1.19 **FAILS**. The text says "the u-p block of the eigenstress and the n_l channel are omitted". At the tag, the live-K chain adds a K[u,p] block of the eigenstress (FEM l.7667–7673, `dsig_sw_dp`), and the level form rescales it by the implicit factor (l.7641–7652). An implicit n_l/ρ_lR channel is added at KKT-active points under `micro_ceiling_sw_tangent = kkt_active` (l.7274–7282). An n_S channel is added (l.7735–7757). The header comment l.7694–7698 still lists the n_S dependence as "NOT included", so it is stale relative to l.7735ff.
*Corrected text:* "The tangent omits the direct n_l (u-p) channel of the eigenstress away from KKT-active points and the u-p entry of the n_S channel; the drain-feedback, n_S-channel and live-K blocks that it does add are labelled DIAGNOSTIC or JACOBIAN-ONLY (l.7635–7673, l.7677–7757)."

1.20 **HOLDS**. CRP l.1130–1146 is one `INFO("DIAGNOSTIC, NOT FOR PRODUCTION (swelling-stress fixes a/b, 2026-10-04): swelling_stress_K_level = {}, swelling_stress_form = {}. ...")`. It covers both switches and is printed only when one of them is on (the enclosing `if`, l.1094–1095).

1.21 **HOLDS**. The clamp guard: FEM l.7460–7472 is the comment `RULED 2026-10-05 ... always on ... reads the phi_prev of the porosity law at EVERY integration point, unconditionally ... JACOBIAN-ONLY: residual untouched.`, and the unconditional read is at l.7473–7477. It sits inside the live-K block (l.7283–7284), so "unconditional" means "at every IP of the live-K chain".

1.22 **HOLDS**. Deviation (1): "in the level form it enters the level itself (l.3572)". l.3572 is `-n_S * n_l * (...)`, and n_S = max(1e-16, 1 − φ_M) (FEM l.3910–3913).

## (2) Statements labelled *derived* (sympy, `audit_level_form_sympy.py`)

2.1 **HOLDS**. σ' = ∂Ψ/∂ε = ∂Ψ^skel/∂ε + ∂Ψ^ads/∂ε follows from additivity. Because Ψ^ads depends on ε only through ε_v, σ_sw is isotropic.

2.2 **CANNOT CHECK**. p_cap = −(1/φ_M) ∂Ψ^cap/∂S_L is a sign and weight convention. It has no code or sympy counterpart in scope.

2.3 **FAILS**. μ_ℓR = ∂Ψ/∂n_ℓ /((1−φ_M)ρ_ℓR) is listed among the conjugates that "follow by the Coleman–Noll argument on ε and S_L". n_ℓ is not one of those arguments, and the appendix itself later calls it an internal variable (A_m = −∂Ψ/∂n_ℓ).
*Corrected text:* "σ' and p_cap follow by the Coleman–Noll argument on ε and S_L; μ_ℓR is defined as the n_ℓ derivative of Ψ divided by (1−φ_M)ρ_ℓR (a definition, since n_ℓ is an internal variable)."

2.4 **HOLDS**. Eq. BGR_dv_film: −∫₀^{ε_v}(1−φ_M)n_ℓ(Π + bKε)dε equals the closed text, and ∂Ψ_film/∂ε_v = −(1−φ_M)n_ℓ(Π + bKε_v) (sympy OK). Note: "drained line" here means σ'_mean = K_sk ε_v, without σ0 and without σ_sw. That is not the σ'_mean of Eq. BGR_app_p, and the appendix acknowledges the difference for the level form.

2.5 **HOLDS**. Eq. BGR_dv_partner: [∂Ψ_film/∂n_ℓ]/((1−φ_M)ρ) = −[(Π + n_ℓΠ')ε_v + ½bKε_v²]/ρ (sympy OK). It also equals PE l.584–587 (sympy OK).

2.6 **HOLDS**. Eq. BGR_dv_ident: ∂σ_sw/∂n_ℓ = (1−φ_M)ρ ∂μ_mech/∂ε_v = −(1−φ_M)(Π + n_ℓΠ' + bKε_v) (sympy OK, both equalities).

2.7 **HOLDS**. Eq. BGR_dv_closed: solving s + L⁽⁰⁾ = −φ_m[Π + b(σ0 + Kε_v + s)] gives s = −[φ_m(Π + b(σ0 + Kε_v)) + L⁽⁰⁾]/(1 + bφ_m), which matches the text (sympy OK). It requires σ_sw⁽⁰⁾ = 0 (see 1.5).

2.8 **HOLDS**. Telescoping makes σ_sw a state function; the algebra is trivial. The unit test `TEST(SwellingStressLevelForm, LevelFormIsAStateFunction)` (Tests/ProcessLib/RichardsMechanics/SwellingStressLevelForm.cpp l.218–255) checks it at Gauss-point level.

2.9 **HOLDS**. Eq. BGR_dv_energy: ∂Ψ^cpl/∂ε_v − L⁽⁰⁾ = σ_sw of Eq. BGR_dv_closed (sympy OK). "less the constant L⁽⁰⁾" is correct.

2.10 **HOLDS**. Eq. BGR_dv_need: [∂Ψ^cpl/∂n_ℓ]/((1−φ_M)ρ) equals the text (sympy OK).

2.11 **HOLDS**. Eq. BGR_dv_need and the branch-tip partner ρ μ_m = −[(Π + b(s0 − L0))ε_v + ½bK_skε_v²]/Z² − n_ℓΠ'ε_v/Z, with Z = 1 + bφ_m, are the same object. The tip's version is in docs/dsm_z_weight_design_2026-10.md l.35 and Tests/Python/dsm_z_weight_partner_sympy.py l.33. The symbols map as n_S ≡ (1−φ_M) and s0 ≡ σ0, and the difference is identically 0 in sympy. Neither one is wrong.

2.12 **HOLDS**. The Maxwell identity for Ψ^cpl, ∂σ_sw/∂n_ℓ = (1−φ_M)ρ ∂μ_need/∂ε_v, closes (sympy OK).

2.13 **HOLDS**. Eq. BGR_dv_miss reproduces exactly (sympy OK), and it vanishes at b = 0. Note: the miss does **not** vanish at ε_v = 0. Its value there is (1−φ_M)[Π + n_ℓΠ' − (Π + b(σ0 − L0))/Z² − n_ℓΠ'/Z].

2.14 **HOLDS**. "Both partners vanish where ε_v = 0" checks out in sympy.

2.15 **CANNOT CHECK**. "...as in the one-element cell of Model I". The path test of the same appendix applies a 1 % axial strain to a Model I cell, so ε_v ≠ 0 on at least one path. The sentence holds only for an isochoric Model I cell, and the cell kinematics are not in the repository.

2.16 **HOLDS**. "The partners differ in the Π terms, in σ0 − L0 (absent in the code) and in the stiffness term" matches the sympy difference: factors 1/Z² on Π and on ½bKε², 1/Z on n_ℓΠ', plus −b(σ0 − L0)ε_v/Z². At the tag, PE l.584–587 has no σ0 − L0 term. The tip's Z-weighted partner (PE l.621) does not exist at the tag.

2.17 **FAILS** (minor). Eq. BGR_dv_dissex: "α_M(Δμ)² ≥ 0 ⟺ α_M ≥ 0" is false for Δμ = 0.
*Corrected text:* "D_exch ≥ 0 for every μ_LR − μ_ℓR if and only if α_M ≥ 0."

2.18 **FAILS**. "The exchange law is the kinetic law that realises it [D = A_m ṅ_ℓ ≥ 0]" is asserted. Two sentences later the appendix says the identification of A_m ṅ_ℓ with (μ_LR − μ_ℓR)ρ̂_ℓ is not derived.
*Corrected text:* "The exchange law is proposed as the kinetic law that realises it; the identification is not derived here."

## (3) Symbol usage

3.1 **HOLDS**. n_S = 1 − φ_M is the weight of the level. FEM l.3910–3913 `n_S = std::max(1e-16, 1.0 - phi_M)` (TransportPorosityData), l.3572, and the comment at l.3598.

3.2 **FAILS**. Eqs. BGR_dv_ident and BGR_dv_miss are taken "at constant (1−φ_M)". The code treats n_S as dependent on n_ℓ at fixed total porosity: n_S = (1 − φ)/(1 − n_ℓ), so dn_S/dn_ℓ = n_S/(1 − n_ℓ) (FEM l.7741–7743, and l.7589–7603, which add L/(1 − n_ℓ) to the n_ℓ partial of the level). The miss of the code's state therefore carries an extra term (1−φ_M)·L/((1−φ_M)(1 − n_ℓ)) that neither partner contains.
*Corrected text:* "The miss is given at constant (1−φ_M); in the code n_S = (1−φ)/(1−n_ℓ) moves with n_ℓ at fixed φ, which adds L/(1−n_ℓ) to ∂σ_sw/∂n_ℓ and is not in either partner."

3.3 **HOLDS** (conditional). In the level, Π is the bare disjoining pressure, Π = −ρ μ_vdW (FEM l.3515–3516, `Pi_curr_film = -p_L_m_density_curr_film * mu_lR_curr_film`). μ is the vdW+augmentation potential without μ_mech (l.3487–3504). The density is ρ_ℓR only when `use_micro_liquid_density_for_micro_pressure = true` (l.3509–3513); otherwise it is ρ_LR (CRP l.862–870). Whether the submitted deck uses ρ_ℓR cannot be read from the repository.

3.4 **FAILS**. The Π inside μ_mech is not always the Π of the level. At the fold, Π = −ρ_ℓR,used·out.mu_lR (FEM l.986–989), and `out` has already been multiplied by the macro-floor cutoff g (FEM l.810, l.760–802; Π' then also contains μ₀ dg/dn_ℓ). The eigenstress Π has no cutoff. The two agree only where `macro_porosity_floor` ≤ 0 or g = 1. Whether the submitted deck sets a floor cannot be checked here.
*Corrected text:* "Π in μ_mech is the cutoff-weighted g·Π (macro-floor cutoff), which equals the Π of the level only where the cutoff does not act."

3.5 **HOLDS**. Sign convention: tension is positive, p_conf = −tr σ'/3 (FEM l.3970–3975), and L = −n_S n_ℓ(Π − b p_conf) = −n_S n_ℓ(Π + b σ'_mean). The swelling stress enters as ε_m = ε + C_el⁻¹σ_sw (FEM l.7622–7624), so σ'_mean = σ0 + Kε_v + σ_sw. Every formula of the appendix carries these signs consistently (sympy chain 2.4–2.13).

3.6 **HOLDS**. ρ_d = ρ_s(1 − φ): the code uses `micro_solid_density_reference * (1.0 - phi)` (PEP l.978–981).

## (4) Measured claims

4.1 **CANNOT CHECK**. 0.46 Pa (expansion, level form). The repository holds no one-element path-test output, and a grep of the whole tag finds none. The numbers would sit in the run directories of the path test on Vinay's machines (e.g. `~/ogs-models/...`), not in the repo.

4.2 **CANNOT CHECK**. 0.74 Pa (compression, level form). Same as 4.1.

4.3 **CANNOT CHECK**. 34.0 kPa (step rule, expansion). Same as 4.1.

4.4 **CANNOT CHECK**. 12.3 kPa (step rule, compression). Same as 4.1.

4.5 **CANNOT CHECK**. The one-element closed loops of 2026-10-13 are not in this repository.

## (5) The last paragraph and the consistency claims

5.1 **HOLDS**. Deviations (2), (4) and (5) stand open at the tag. μ_mech (PE l.584–587) has no drain (σ0 − L0 or b σ'_mean) term and no dK/dρ_d term, and it is the drained-line form.

5.2 **HOLDS**. "1 stands open" is consistent with the code: n_S moves with n_ℓ (3.2) and no partner carries that dependence. The appendix gives a weaker reason, though: "it enters the level itself".

5.3 **HOLDS**. Deviation (3), the path of K_aug, is closed for the stress by the per-level K (FEM l.3258–3262, PEP l.781–785). The identities remain at constant K_aug, which the appendix states as open.

5.4 **HOLDS**. Item (i) of the open list: below the micro maximum, the partner enters every micro solve and the exchange through the fold at FEM l.1002. At KKT-active points the booked exchange is S_s (l.1998), so the partner does not enter the booked exchange there.

5.5 **HOLDS**. Item (ii): at the tag there is no code for μ_need. `computeZWeightedMechanicalMicroPotential` exists only at the branch tip (PE l.621 at cc5e71e12a).

5.6 **HOLDS**. Items (iii) and (iv) are consistent with the code. Only the drained-line partner is implemented, and Ψ^cpl requires a linear skeleton (the correction loop at FEM l.3784ff. is for non-linear skeletons).

5.7 **FAILS**. The appendix implies thermodynamic consistency of the submitted route in one sentence: "The film reservoir carries ε_v, so that the swelling eigenstress is its strain conjugate and is not imposed." This states as fact, for the eigenstress of the reported runs, what holds only on the drained line with the film frozen. The appendix's own Eq. BGR_dv_miss shows that for the level form the conjugacy is missed. Elsewhere the text is scoped correctly ("open", "not shown", "drained line").
*Corrected text:* "The film reservoir carries ε_v, so that on the drained line with the film frozen the swelling eigenstress is its strain conjugate; for the level form of the submitted runs it is not (Eq. BGR_dv_miss)."

5.8 **FAILS**. "The exchange law is the kinetic law that realises it" implies a dissipation inequality for the submitted exchange that is not derived (see 2.18). At KKT points the booked exchange is not the law (FEM l.1987–1998); the appendix states that correctly elsewhere.
*Corrected text:* as 2.18.

## Count

HOLDS 43 · FAILS 8 · CANNOT CHECK 8 (total 59).

FAILS: 1.19, 2.3, 2.17, 2.18, 3.2, 3.4, 5.7, 5.8.
