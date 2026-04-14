BeginPackage["DsmMicromacroDriver`"];

CreateDefaultParameters::usage =
  "CreateDefaultParameters[] returns an Association of default DSM driver parameters.";

CreateDefaultInitialState::usage =
  "CreateDefaultInitialState[params] returns a physically bounded initial microstate.";

RunHistory::usage =
  "RunHistory[steps, params, initialState] advances the local DSM micro-macro state over a step list.";

RunDryDensitySweep::usage =
  "RunDryDensitySweep[dryDensityList, params, opts] runs repeated local-driver histories over dry-density cases.";

ExportHistory::usage =
  "ExportHistory[result, basePath] writes result rows to basePath.csv and full data to basePath.json.";

Begin["`Private`"];

toNumeric[value_, fallback_] := If[NumericQ[value], N[value], fallback];
clip[value_, bounds_List] := Clip[N[value], bounds];

CreateDefaultParameters[] :=
  <|
    "PressureTolerance" -> 1.*^-12,
    "PorosityUpper" -> 1. - 1.*^-12,
    "NFloor" -> 1.*^-16,
    "RhoFloor" -> 1.*^-12,
    "RhoLRMacro" -> 1000.,
    "RhoL0Micro" -> 1300.,
    "RhoSR" -> 2470.,
    "DensityA" -> 1.3,
    "DensityB" -> 1.,
    "HamakerConstant" -> -6.*^-20,
    "SpecificSurface" -> 100.,
    "MicroPotentialSign" -> 1.,
    "VdwMultiplier" -> 1.,
    "MassExchangeCoefficient" -> 1.,
    "MacroViscosity" -> 1.,
    "SwellingSlope" -> 0.1,
    "BulkModulus" -> 6.666666666666667*^9,
    "UseKinematicPorosityUpdate" -> True,
    "AreaFactorTuller" -> 1.,
    "PoreAreaShapeFactorTuller" -> 0.8584073464102069,
    "CharacteristicPoreSize" -> 1.*^-5,
    "SurfaceTension" -> 0.0715,
    "MechanicalMode" -> "MCCCarrier"
  |>;

CreateDefaultInitialState[params_: Automatic] := Module[
  {p, phiTotal, phiM, nL0, rho0},
  p = If[AssociationQ[params], params, CreateDefaultParameters[]];
  phiTotal = 0.2;
  nL0 = 0.059382805374452955;
  phiM = Max[0., phiTotal - nL0];
  rho0 = 2193.3876099104677;
  <|
    "n_l" -> clip[nL0, {p["NFloor"], phiTotal}],
    "rho_lR" -> Max[p["RhoFloor"], rho0],
    "phi_m" -> clip[nL0, {p["NFloor"], phiTotal}],
    "phi_M" -> phiM,
    "epsilon_sw" -> -0.004061719462554706,
    "sigma_S_xx" -> 2.7078129750364702*^7,
    "sigma_S_yy" -> 2.7078129750364702*^7,
    "sigma_S_zz" -> 2.7078129750364702*^7
  |>
];

saturationFromTuller[pressure_, params_Association] := Module[
  {ptol, areaFactor, shapeFactor, poreSize, gamma, capillaryTerm},
  ptol = params["PressureTolerance"];
  areaFactor = params["AreaFactorTuller"];
  shapeFactor = params["PoreAreaShapeFactorTuller"];
  poreSize = params["CharacteristicPoreSize"];
  gamma = params["SurfaceTension"];
  If[pressure > -ptol,
    1.,
    capillaryTerm =
      (4. * shapeFactor * gamma^2) / (areaFactor * poreSize^2 * pressure^2);
    1. - Exp[-capillaryTerm]
  ]
];

muMacroFromPressure[pressure_, rhoLR_, ptol_] :=
  If[pressure > -ptol, 0., pressure / rhoLR];

omegaFromState[nL_, rhoLr_, nS_, rhoSR_] :=
  Max[1.*^-30, nL * rhoLr / Max[1.*^-30, nS * rhoSR]];

rhoMicroEOS[nL_, rhoLrGuess_, params_Association, nS_] := Module[
  {rhoL0, rhoLR, aRho, bRho, rhoSR, omega},
  rhoL0 = params["RhoL0Micro"];
  rhoLR = params["RhoLRMacro"];
  aRho = params["DensityA"];
  bRho = params["DensityB"];
  rhoSR = params["RhoSR"];
  omega = omegaFromState[nL, rhoLrGuess, nS, rhoSR];
  rhoL0 * Exp[-aRho * omega^bRho] + rhoLR
];

muMicroFromState[nL_, rhoLr_, params_Association, nS_, vdwMultiplier_: 1.] := Module[
  {hamakerMagnitude, signFactor, specificSurface, rhoSR, omega},
  hamakerMagnitude = Abs[params["HamakerConstant"]];
  signFactor = params["MicroPotentialSign"] * Sign[params["HamakerConstant"]];
  specificSurface = params["SpecificSurface"];
  rhoSR = params["RhoSR"];
  omega = omegaFromState[nL, rhoLr, nS, rhoSR];
  signFactor * vdwMultiplier * hamakerMagnitude * rhoLr^3 / (6. * Pi) *
    (specificSurface / omega)^3
];

computeTotalPorosityTrial[state_Association, deltaEpsV_, params_Association] := Module[
  {phiPrev, phiTrial, denom},
  phiPrev = clip[
    toNumeric[state["phi_m"], 0.] + toNumeric[state["phi_M"], 0.],
    {params["NFloor"], params["PorosityUpper"]}
  ];
  If[TrueQ[params["UseKinematicPorosityUpdate"]],
    denom = 1. + deltaEpsV;
    phiTrial =
      If[NumericQ[denom] && Abs[denom] > 1.*^-12,
        (phiPrev + deltaEpsV) / denom,
        phiPrev
      ],
    phiTrial = phiPrev
  ];
  clip[phiTrial, {params["NFloor"], params["PorosityUpper"]}]
];

solveRhoForNl[nL_, rhoInit_, params_Association, nS_] := Module[
  {rho = Max[params["RhoFloor"], rhoInit], rhoNext, iter},
  For[iter = 1, iter <= 100, iter++,
    rhoNext = Max[params["RhoFloor"], rhoMicroEOS[nL, rho, params, nS]];
    If[Abs[rhoNext - rho] <= 1.*^-12 * Max[1., Abs[rho]], Return[rhoNext]];
    rho = rhoNext;
  ];
  rho
];

solveMicrostateStep[state_Association, step_Association, params_Association] := Module[
  {
    dt, pressure, deltaEpsV, phiTotalTrial, nUpper, nFloor, rhoFloor, nS,
    nPrev, rhoPrev, prevMass, volRate, rhoLR, muLR, alphaEff, vdwMultiplier,
    root, nVar, rhoVar, nCandidate, rhoCandidate, converged, residual1,
    residual2, muMicro, rhoHat, nExplicit, rhoResolved
  },
  dt = Max[0., toNumeric[Lookup[step, "dt", 1.], 1.]];
  pressure = toNumeric[Lookup[step, "pressure", 0.], 0.];
  deltaEpsV = toNumeric[Lookup[step, "delta_epsilon_v", 0.], 0.];
  vdwMultiplier = toNumeric[Lookup[step, "vdw_multiplier", params["VdwMultiplier"]], params["VdwMultiplier"]];

  nFloor = params["NFloor"];
  rhoFloor = params["RhoFloor"];
  phiTotalTrial = computeTotalPorosityTrial[state, deltaEpsV, params];
  nUpper = Max[nFloor, phiTotalTrial];
  nS = Max[nFloor, 1. - phiTotalTrial];

  nPrev = clip[toNumeric[state["n_l"], nFloor], {nFloor, nUpper}];
  rhoPrev = Max[rhoFloor, toNumeric[state["rho_lR"], params["RhoLRMacro"]]];
  prevMass = nPrev * rhoPrev;
  volRate = If[dt > 0., deltaEpsV / dt, 0.];
  rhoLR = params["RhoLRMacro"];
  muLR = muMacroFromPressure[pressure, rhoLR, params["PressureTolerance"]];
  alphaEff = params["MassExchangeCoefficient"] * rhoLR / params["MacroViscosity"];

  converged = False;
  root = Quiet @ Check[
    FindRoot[
      {
        nVar * rhoVar - prevMass -
          dt * alphaEff * (muLR - muMicroFromState[nVar, rhoVar, params, nS, vdwMultiplier]) -
          dt * nVar * rhoVar * volRate == 0,
        rhoVar - rhoMicroEOS[nVar, rhoVar, params, nS] == 0
      },
      {
        {nVar, SetPrecision[nPrev, 50]},
        {rhoVar, SetPrecision[rhoPrev, 50]}
      },
      WorkingPrecision -> 50,
      AccuracyGoal -> 20,
      PrecisionGoal -> 20,
      MaxIterations -> 100
    ],
    $Failed,
    {FindRoot::lstol, FindRoot::cvmit, FindRoot::jsing, FindRoot::nlnum}
  ];

  If[root =!= $Failed && NumericQ[nVar /. root] && NumericQ[rhoVar /. root],
    nCandidate = clip[N[nVar /. root], {nFloor, nUpper}];
    rhoCandidate = Max[rhoFloor, N[rhoVar /. root]];
    residual1 =
      nCandidate * rhoCandidate - prevMass -
        dt * alphaEff * (muLR - muMicroFromState[nCandidate, rhoCandidate, params, nS, vdwMultiplier]) -
        dt * nCandidate * rhoCandidate * volRate;
    residual2 = rhoCandidate - rhoMicroEOS[nCandidate, rhoCandidate, params, nS];
    converged = Abs[residual1] < 1.*^-8 && Abs[residual2] < 1.*^-8;
  ];

  If[!converged,
    muMicro = muMicroFromState[nPrev, rhoPrev, params, nS, vdwMultiplier];
    rhoHat = alphaEff * (muLR - muMicro);
    nExplicit = nPrev + dt * rhoHat / Max[rhoFloor, rhoPrev] + dt * nPrev * volRate;
    nCandidate = clip[nExplicit, {nFloor, nUpper}];
    rhoResolved = solveRhoForNl[nCandidate, rhoPrev, params, nS];
    rhoCandidate = Max[rhoFloor, rhoResolved];
  ];

  muMicro = muMicroFromState[nCandidate, rhoCandidate, params, nS, vdwMultiplier];
  rhoHat = alphaEff * (muLR - muMicro);

  <|
    "n_l" -> nCandidate,
    "rho_lR" -> rhoCandidate,
    "mu_LR" -> muLR,
    "mu_lR" -> muMicro,
    "rho_l_hat" -> rhoHat,
    "phi_total_trial" -> phiTotalTrial,
    "n_upper_bound" -> nUpper,
    "solver_status" -> If[converged, "newton", "fallback_explicit"]
  |>
];

advanceState[state_Association, step_Association, params_Association, stepIndex_Integer] := Module[
  {
    deltaEpsV, pressure, dt, micro, phiM, phiMPrev, epsSwPrev,
    deltaEpsSw, epsSw, sigmaPrev, sigmaNew, satMacro, phiTotal, nS, nL
  },
  deltaEpsV = toNumeric[Lookup[step, "delta_epsilon_v", 0.], 0.];
  pressure = toNumeric[Lookup[step, "pressure", 0.], 0.];
  dt = Max[0., toNumeric[Lookup[step, "dt", 1.], 1.]];
  phiMPrev = toNumeric[state["phi_m"], 0.];
  epsSwPrev = toNumeric[state["epsilon_sw"], 0.];
  sigmaPrev = toNumeric[state["sigma_S_xx"], 0.];

  micro = solveMicrostateStep[state, step, params];
  phiM = Max[0., micro["phi_total_trial"] - micro["n_l"]];
  phiTotal = micro["n_l"] + phiM;
  nS = Max[params["NFloor"], 1. - phiTotal];
  satMacro = saturationFromTuller[pressure, params];
  nL = phiM * satMacro;

  deltaEpsSw = params["SwellingSlope"] * (micro["n_l"] - phiMPrev);
  epsSw = epsSwPrev + deltaEpsSw;
  sigmaNew =
    sigmaPrev + params["BulkModulus"] * deltaEpsV -
      params["BulkModulus"] * deltaEpsSw;

  <|
    "step" -> stepIndex,
    "pressure" -> pressure,
    "dt" -> dt,
    "delta_epsilon_v" -> deltaEpsV,
    "S_L" -> satMacro,
    "mu_LR" -> micro["mu_LR"],
    "n_l" -> micro["n_l"],
    "phi_m" -> micro["n_l"],
    "phi_M" -> phiM,
    "phi" -> phiTotal,
    "n_S" -> nS,
    "n_L" -> nL,
    "rho_lR" -> micro["rho_lR"],
    "rho_LR" -> params["RhoLRMacro"],
    "omega_l" -> omegaFromState[micro["n_l"], micro["rho_lR"], nS, params["RhoSR"]],
    "mu_lR" -> micro["mu_lR"],
    "rho_l_hat" -> micro["rho_l_hat"],
    "delta_epsilon_sw" -> deltaEpsSw,
    "epsilon_sw" -> epsSw,
    "sigma_S_xx" -> sigmaNew,
    "sigma_S_yy" -> sigmaNew,
    "sigma_S_zz" -> sigmaNew,
    "sigma_S_xy" -> 0.,
    "sigma_S_yz" -> 0.,
    "sigma_S_xz" -> 0.,
    "solver_status" -> micro["solver_status"],
    "n_upper_bound" -> micro["n_upper_bound"],
    "mechanical_mode" -> params["MechanicalMode"]
  |>
];

RunHistory[steps_List, params_: Automatic, initialState_: Automatic] := Module[
  {p, state, rows = {}, i, stepAssoc, row, stateKeys},
  p = If[AssociationQ[params], params, CreateDefaultParameters[]];
  state = If[AssociationQ[initialState], initialState, CreateDefaultInitialState[p]];
  stateKeys = {"n_l", "rho_lR", "phi_m", "phi_M", "epsilon_sw", "sigma_S_xx", "sigma_S_yy", "sigma_S_zz"};

  For[i = 1, i <= Length[steps], i++,
    stepAssoc = If[AssociationQ[steps[[i]]], steps[[i]], <|"pressure" -> toNumeric[steps[[i]], 0.]|>];
    row = advanceState[state, stepAssoc, p, i - 1];
    AppendTo[rows, row];
    state = Join[state, KeyTake[row, stateKeys]];
  ];

  <|
    "rows" -> rows,
    "final_state" -> state,
    "parameters" -> p
  |>
];

resolveMultiplier[dryDensity_, params_Association, mode_] := Module[
  {rules, xVals, yVals, interp},
  Which[
    NumericQ[mode], N[mode],
    Head[mode] === Function, N[mode[dryDensity]],
    AssociationQ[mode] && Length[mode] > 1,
      rules = SortBy[Normal[mode], First];
      xVals = rules[[All, 1]];
      yVals = rules[[All, 2]];
      interp = Interpolation[Transpose[{xVals, yVals}], InterpolationOrder -> 1];
      N[interp[dryDensity]],
    AssociationQ[mode] && KeyExistsQ[mode, dryDensity], N[mode[dryDensity]],
    True, N[params["VdwMultiplier"]]
  ]
];

Options[RunDryDensitySweep] = {
  "SolidDensity" -> 2700.,
  "InitialMicroWaterFraction" -> 0.99,
  "CaseSteps" -> 120,
  "Pressure" -> 0.,
  "DeltaEpsilonV" -> 0.,
  "Dt" -> 1.,
  "MultiplierByDryDensity" -> Automatic
};

RunDryDensitySweep[dryDensityList_List, params_: Automatic, opts : OptionsPattern[]] := Module[
  {
    p, solidDensity, microWaterFraction, caseSteps, pressure, deltaEpsV, dt,
    multiplierMode, rows = {}, dd, phiTotal, phiM0, init, caseParams, steps,
    result, finalRow, solverFallbackCount, multiplier
  },
  p = If[AssociationQ[params], params, CreateDefaultParameters[]];
  solidDensity = OptionValue["SolidDensity"];
  microWaterFraction = OptionValue["InitialMicroWaterFraction"];
  caseSteps = Round[OptionValue["CaseSteps"]];
  pressure = OptionValue["Pressure"];
  deltaEpsV = OptionValue["DeltaEpsilonV"];
  dt = OptionValue["Dt"];
  multiplierMode = OptionValue["MultiplierByDryDensity"];

  Do[
    phiTotal = clip[1. - dd / solidDensity, {p["NFloor"], p["PorosityUpper"]}];
    phiM0 = clip[microWaterFraction * phiTotal, {p["NFloor"], phiTotal}];
    multiplier = resolveMultiplier[dd, p, multiplierMode];

    caseParams = Join[p, <|"VdwMultiplier" -> multiplier|>];
    init = <|
      "n_l" -> phiM0,
      "rho_lR" -> Max[p["RhoFloor"], p["RhoL0Micro"]],
      "phi_m" -> phiM0,
      "phi_M" -> Max[0., phiTotal - phiM0],
      "epsilon_sw" -> 0.,
      "sigma_S_xx" -> 0.,
      "sigma_S_yy" -> 0.,
      "sigma_S_zz" -> 0.
    |>;
    steps = Table[
      <|
        "pressure" -> pressure,
        "delta_epsilon_v" -> deltaEpsV,
        "dt" -> dt,
        "vdw_multiplier" -> multiplier
      |>,
      {caseSteps}
    ];

    result = RunHistory[steps, caseParams, init];
    finalRow = Last[result["rows"]];
    solverFallbackCount =
      Count[Lookup[result["rows"], "solver_status", ""], "fallback_explicit"];

    AppendTo[
      rows,
      <|
        "dry_density_kg_m3" -> dd,
        "phi_total_initial" -> phiTotal,
        "vdw_multiplier" -> multiplier,
        "final_n_l" -> finalRow["n_l"],
        "final_phi_m" -> finalRow["phi_m"],
        "final_phi_M" -> finalRow["phi_M"],
        "final_phi" -> finalRow["phi"],
        "final_rho_lR" -> finalRow["rho_lR"],
        "final_mu_lR" -> finalRow["mu_lR"],
        "final_rho_l_hat" -> finalRow["rho_l_hat"],
        "final_swelling_stress_mpa" -> finalRow["sigma_S_xx"] / 1.*^6,
        "fallback_step_count" -> solverFallbackCount
      |>
    ];
    ,
    {dd, dryDensityList}
  ];

  <|
    "rows" -> rows,
    "parameters" -> p,
    "options" -> <|
      "SolidDensity" -> solidDensity,
      "InitialMicroWaterFraction" -> microWaterFraction,
      "CaseSteps" -> caseSteps,
      "Pressure" -> pressure,
      "DeltaEpsilonV" -> deltaEpsV,
      "Dt" -> dt
    |>
  |>
];

ExportHistory[result_Association, basePath_String] := Module[
  {rows, csvPath, jsonPath, header},
  rows = Lookup[result, "rows", {}];
  csvPath = basePath <> ".csv";
  jsonPath = basePath <> ".json";
  If[rows === {},
    Export[csvPath, {{}}, "CSV"],
    header = Keys[First[rows]];
    Export[
      csvPath,
      Prepend[(Lookup[#, header] &) /@ rows, header],
      "CSV"
    ]
  ];
  Export[jsonPath, result, "RawJSON"];
  <|"csv" -> csvPath, "json" -> jsonPath|>
];

End[];

EndPackage[];
