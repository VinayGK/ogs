# ANCHORS MS33 — FASTEX recalibration provenance (2026-06-02)

Working folder for the exchange-rate + force-augmentation recalibration with a
constant macro relative permeability. **Temporary** (constant k_rel) — reverts to an
orthotropic power-law k_rel later (Vinay).

## Motivation (verified)

Model III/IV centre stresses did **not** flatten within the 200 d window (creep past
200 d), unlike the inter-team results. Diagnosis: the macro stays dry (S_L≈0 in the
bulk; only the bottom wets), so with the original k_rel=S_e³ floored at 1e-2 the
**dry-macro relative permeability collapsed (k_r≈0.01)** and throttled water transport
up the column — the micro filled as a slow front and the micro pressure relaxed past
200 d. The limiter was hydraulic transport, not the exchange rate.

Separately, the Model I "flat" plateaus were **not** true equilibrium: identical micro
water content (0.4964) but α-dependent **micro pressure** (12.0 vs 15.87 MPa for
α=1e-13 vs 3e-13). The original K was therefore anchored to a frozen, sub-equilibrium
micro pressure.

## Calibration design (§2-clean: two observables, two knobs)

- **k_rel (constant 0.1)** — fixes the dry-macro transport bottleneck (≥10% conductivity
  everywhere). Applied to ALL models. *Temporary simplification.*
- **α_exchange = `mass_exchange_coefficient` = 3e-13** — calibrated to the **transient**:
  Model I t95 = 18/12/10 d (dd1400/1600/1800), all <20 d, and III flattens (~50 d, flat
  to 200 d) instead of creeping.
- **K = `vdw_augmentation_prefactor`** — calibrated to the **Dixon EMDD≡ρ_d magnitude**
  at the fixed α (single linear scale, verified):

| dd | K (J/kg) | p_ss (MPa) | Dixon target | err |
|---|---|---|---|---|
| 1400 | 34368 | 4.917 | 4.922 | −0.1% |
| 1600 | 83377 | 14.166 | 14.161 | +0.0% |
| 1800 | 224610 | 40.859 | 40.860 | −0.0% |

III/IV/VII inherit the dd1.6 regime: K scaled by 0.97729 (85312.6→83377), preserving
IV's clay/pellet ratio.

Order matters: α fixes *where* the micro pressure settles, so α is set first (transient),
then K is fit *at that α* (magnitude). No chasing α→∞; K is calibrated at the production α.

## Calibration anchor (CLAUDE.md §12.1)

Dixon et al. (2023), Applied Clay Science 241:106998, MX-80 Fig. 1 median,
σ_swell[MPa] = 0.003·exp(5.2883·EMDD[Mg/m³]), EMDD ≡ ρ_d (working-group agreement
2026-05-27). Targets 4.922 / 14.161 / 40.860 MPa at ρ_d = 1.4 / 1.6 / 1.8 g/cm³.

## Known limitations (documented, not hidden)

- **dd1800 saturated-corner singularity**: with constant k_rel=0.1 the saturated macro
  conductivity drops 10×, ill-conditioning the high-K hydro-mechanical coupling at the
  20 d ramp→constant corner. OGS aborts in Eigen SparseLU initialization (singular
  matrix) even at dt=1e-4 — **structural, not a dt floor** (dt-refine `minimum_dt`→1e-4
  did not help). The 0–20 d path + the Dixon endpoint (40.86 MPa, flat by 10 d) are
  captured. **RESOLVED (agent rerun 2026-06-02):** the power-law k_rel (min 0.1) on the
  single element completes the full 200 d at K=224610 -> 40.859 MPa (Dixon, -0.0%), **no K
  re-fit needed** (the power law restores k_r=1 at saturation, so the K-set endpoint is
  unchanged). See `LE/I_dd1800_pl/` + `run_fastex_pl.py`. Model I needs no transport floor
  (single element), so this power-law exception is clean.
- **MCC suite (agent rerun confirmed):** the k_rel fix made **MCC III + IV converge to
  200 d** (they crashed at the apex in the prior deliverable) — and **MCC IV homogenises**
  (clay 1585->1408, pellet 900->1016, BExM-style; the LE skeleton does not). MCC dd1600
  **yields** to ~11.8-14.8 MPa (genuine cam-clay plasticity at dd1.6; non-uniform field).
  MCC dd1800 (return-map at K=224610, status -1, t~0.04 d) and MCC VII (tension apex,
  p_net->0 with residual q, status -1, t=10 s) are **genuinely upstream MFront** return-map
  failures, NOT k_rel — the power-law variant did not help them. MCC Model I dd1400 ≡ LE.

## Files

`run_fastex.py` (prepare+run one model: krel=0.1, α, K set/scale, min_dt=1e-4, meshes
copied local), `run_suite.py` (LE+MCC, 4-parallel), `calib_fastex.py` (the α probe),
`SUITE_SUMMARY.txt` (results, written on suite completion).
