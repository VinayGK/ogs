ClearAll[
  here, outputDir, calibrationCsvPath, params, initialState,
  overlapSteps, overlapResult, overlapExport,
  strainSteps, strainResult, strainExport,
  calibrationRaw, calibrationHeader, ddCol, multCol, multCurve,
  ddList, sweepResult, sweepExport, summaryPath, summary
];

here = DirectoryName[$InputFileName];
outputDir = FileNameJoin[{here, "_outputs"}];
If[!DirectoryQ[outputDir], CreateDirectory[outputDir]];

Get[FileNameJoin[{here, "DsmMicromacroDriver.wl"}]];

params = DsmMicromacroDriver`CreateDefaultParameters[];
initialState = DsmMicromacroDriver`CreateDefaultInitialState[params];

overlapSteps = Table[
  <|
    "pressure" -> 250. * i,
    "delta_epsilon_v" -> 0.,
    "dt" -> 1.
  |>,
  {i, 0, 4}
];

overlapResult =
  DsmMicromacroDriver`RunHistory[overlapSteps, params, initialState];
overlapExport =
  DsmMicromacroDriver`ExportHistory[
    overlapResult,
    FileNameJoin[{outputDir, "overlap_history"}]
  ];

strainSteps = Table[
  <|
    "pressure" -> 1000.,
    "delta_epsilon_v" -> 1.*^-4,
    "dt" -> 1.
  |>,
  {5}
];

strainResult =
  DsmMicromacroDriver`RunHistory[strainSteps, params, overlapResult["final_state"]];
strainExport =
  DsmMicromacroDriver`ExportHistory[
    strainResult,
    FileNameJoin[{outputDir, "strain_coupled_history"}]
  ];

calibrationCsvPath = FileNameJoin[
  {here, "..", "ANCHORS_MS33_ModelI", "villar_dense_dd_calibration.csv"}
];

If[FileExistsQ[calibrationCsvPath],
  calibrationRaw = Import[calibrationCsvPath, "CSV"];
  calibrationHeader = First[calibrationRaw];
  ddCol = First @ FirstPosition[calibrationHeader, "dry_density_kg_m3"];
  multCol = First @ FirstPosition[calibrationHeader, "vdw_multiplier"];
  multCurve =
    Association @ Map[
      (N[ToExpression[#[[ddCol]]]] -> N[ToExpression[#[[multCol]]]]) &,
      Rest[calibrationRaw]
    ],
  multCurve = <|
    1400. -> 1.*^5,
    1450. -> 3.*^5,
    1500. -> 9.*^5,
    1600. -> 8.*^6,
    1700. -> 8.*^7,
    1800. -> 9.*^8
  |>
];

ddList = Sort[Keys[multCurve]];
sweepResult =
  DsmMicromacroDriver`RunDryDensitySweep[
    ddList,
    params,
    "SolidDensity" -> 2700.,
    "InitialMicroWaterFraction" -> 0.99,
    "CaseSteps" -> 120,
    "Pressure" -> 0.,
    "DeltaEpsilonV" -> 0.,
    "Dt" -> 1.,
    "MultiplierByDryDensity" -> multCurve
  ];

sweepExport =
  DsmMicromacroDriver`ExportHistory[
    sweepResult,
    FileNameJoin[{outputDir, "dry_density_sweep"}]
  ];

summaryPath = FileNameJoin[{outputDir, "run_summary.json"}];
summary = <|
  "project" -> "Mathematica_DSM_Driver",
  "overlap_rows" -> Length[overlapResult["rows"]],
  "strain_rows" -> Length[strainResult["rows"]],
  "dry_density_cases" -> Length[sweepResult["rows"]],
  "calibration_curve_source" ->
    If[FileExistsQ[calibrationCsvPath], calibrationCsvPath, "internal_fallback_curve"],
  "exports" -> <|
    "overlap" -> overlapExport,
    "strain" -> strainExport,
    "dry_density_sweep" -> sweepExport
  |>
|>;

Export[summaryPath, summary, "RawJSON"];

Print["Mathematica DSM driver demo completed."];
Print["Overlap history CSV: ", overlapExport["csv"]];
Print["Strain history CSV: ", strainExport["csv"]];
Print["Dry-density sweep CSV: ", sweepExport["csv"]];
Print["Summary JSON: ", summaryPath];

Exit[0];
