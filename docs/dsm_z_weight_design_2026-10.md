# The Z-weighted strain partner of the level-form swelling stress (design note, 2026-10-13)

Status words: posited (statement about the material), derived (follows by the step shown, by symbol),
implemented (code at the named site), computed (a run or a script), open. No claim of thermodynamic
consistency is made; the scope is stated at the end. Branch `dsm_z_weight_2026-10`, cut from tag
`submission-config-anchors-T4.1-T4.2-2026-10-09` (commit 8bbcace35c). Brief:
`ANCHORS/runs/verification/2026-10-13_closed_loop_single_element/BRIEF_Z_weight_cloud_2026-10-13.md`
(ruled with Vinay 2026-10-13 15:55: ONE partner, the level-form partner with divisor Z, s_0 kept).

Symbols. n_l micro water content; phi_M macro porosity; n_S = 1 - phi_M (REV solid fraction, the
level's weight); phi_m = n_S n_l; Pi(n_l) disjoining pressure of the augmented law, Pi' = dPi/dn_l;
b Biot coefficient; K_sk drained bulk modulus of the skeleton; eps_v total volumetric strain; s_0 initial
mean effective stress; L_0 level of the initial state; rho_lR micro water density. OGS sign (tension
positive) throughout; Pi > 0 is repulsive, sigma_sw < 0 is compressive.

## 1. The stress as implemented (implemented)
Level form, `RichardsMechanicsFEM-impl.h`, block `if (level_form_sw)` of
`computeReferenceMicroPorositySwellingStressIncrement` (search `L (1 + c)`):
    L = -n_S n_l [Pi + b sigma'_mean],   sigma_sw = L - L_0,
and, with the elastic prediction m_hat = s_0 + K_sk eps_v of the mean effective stress without the
swelling stress, the closed form  L (1 + c) = -n_S n_l Pi - c (m_hat - L_0),  c = n_S n_l b.

## 2. Solved for the stress (derived)
    sigma_sw (1 + b phi_m) = -phi_m [Pi + b (s_0 - L_0) + b K_sk eps_v]  -  L_0 (1 + b phi_m)
    sigma_sw = -(phi_m/Z) [Pi + b (s_0 - L_0) + b K_sk eps_v] - L_0,     Z = 1 + b phi_m.
The divisor Z is the feedback of the swelling stress on the stress the film feels (gain b phi_m).
Checked against the FEM function in `DSMZWeightedPartnerMaxwellIdentity` (sanity block).

## 3. The energy and the partner (derived)
Assumptions: linear skeleton with constant K_sk; n_S, rho_lR, K_aug, s_0, L_0 held fixed during the
integration; Pi = Pi(n_l) as the code evaluates it (floors and cutoff included through the fold).
Integration over eps_v at fixed n_l, Psi(0) = 0:
    Psi = -(phi_m/Z) [(Pi + b (s_0 - L_0)) eps_v + 1/2 b K_sk eps_v^2]          (dPsi/deps_v = sigma_sw + L_0)
Partner rho_lR mu_m = (1/n_S) dPsi/dn_l, with d(phi_m/Z)/dn_l = n_S/Z^2:
    rho_lR mu_m = -[(Pi + b (s_0 - L_0)) eps_v + 1/2 b K_sk eps_v^2] / Z^2  -  n_l Pi' eps_v / Z.
Maxwell identity d sigma_sw/d n_l = n_S rho_lR d mu_m/d eps_v: closes at every state (sympy
`Tests/Python/dsm_z_weight_partner_sympy.py`, checks 1 to 4 PASS; Google test above, FD against analytic).
At b = 0 (Z = 1) and s_0 = L_0 the partner is the drained-line partner of the tag,
    rho_lR mu_code = -[(Pi + n_l Pi') eps_v + 1/2 b K_sk eps_v^2].

## 4. Difference to the partner of the tag (derived; sign from sympy check 5)
    rho_lR (mu_Z - mu_code) = Pi eps_v (1 - 1/Z^2) + n_l Pi' eps_v (1 - 1/Z) + 1/2 b K_sk eps_v^2 (1 - 1/Z^2)
                              - b (s_0 - L_0) eps_v / Z^2.
For Pi eps_v > 0 the Z partner is the smaller in magnitude (the drained-line partner over-states the
strain response by the factor 1/(1/Z^2) on the Pi term: at 1600 kg/m3 on the micro maximum, phi_m = 0.424,
b = 1, 1/Z^2 = 0.49). Both partners vanish at eps_v = 0.

## 5. Implementation (implemented)
- `PotentialExchangeParameters.h`: `bool micro_potential_partner_z_weight = false` (PRJ tag of the same
  name under `<potential_exchange>`; default off = bit-identical to the tag, no code path changed when off).
- `CreateRichardsMechanicsProcess.cpp`: parse; refused unless `swelling_stress_form = level`, operational
  film route, `film_strain_coupling = off`; INFO line "DIAGNOSTIC, NOT ADOPTED".
- `ConstitutiveRelations/PotentialExchange.h`: `computeZWeightedMechanicalMicroPotential` (value and the
  three tangents d/deps_v, d/dn_l total with n_S held, d/drho_lR) and the dispatcher
  `computeMechanicalMicroPotentialPartner(z_weight, ...)`.
- `RichardsMechanicsFEM-impl.h`: the context carries `rev_solid_fraction` (n_S of the current state, floored
  1e-16 as the swelling correction loop) and `level_reference_offset` = s_0 - L_0, read from the previous
  accepted state as mean(sigma'_prev) - K_sk eps_v,prev - L_prev (`levelFormPartnerInputs`; equal to
  s_0 - L_0 for a linear skeleton at every accepted state). The dispatcher replaces the drained-line call at
  the ONE residual fold point (`applyFilmPressureMicroPotential`) and at the three Jacobian sites (macro
  exchange in assemble, in assembleWithJacobian, and the p-u tangent), so the residual and the tangent see
  the same partner.
- Tests: `DSMZWeightedPartnerMaxwellIdentity` (Google test), `Tests/Python/dsm_z_weight_partner_sympy.py`.

## 6. Scope and what stays open (open)
- Linear skeleton, constant K_sk; n_S, rho_lR and K_aug held in the derivation. The live K(rho_d) table and
  the strain dependence of n_S have no partner here (as before).
- s_0 - L_0 is taken from the previous accepted state; inside a Newton iteration it is therefore constant
  over the step (consistent with the closed form's elastic prediction). Not derived for Cam-clay.
- Whether this partner is the physics, and whether the ANCHORS reruns use it, is Vinay's call. The loop
  test of section C of the brief measures its effect on the stress-work residual; see the run card in the
  ANCHORS lane `z_weight_2026-10-13/`.
