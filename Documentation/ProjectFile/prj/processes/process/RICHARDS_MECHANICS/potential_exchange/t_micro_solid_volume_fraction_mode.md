Selects how the active micro-scale solid volume fraction \f$n_S\f$ (the
code's "active_nS") entering the micro-liquid-density EOS and the
disjoining-pressure / van-der-Waals potential is obtained. Parsed by
`parseMicroSolidVolumeFractionMode()`
(`ProcessLib/RichardsMechanics/CreateRichardsMechanicsProcess.cpp`); an
unrecognised string
`OGS_FATAL("RichardsMechanics: unsupported potential_exchange
micro_solid_volume_fraction_mode '{}'. Currently supported: 'reference',
'current_porosity_split'.")`.

Type: `std::string`, converted to the enum
`ProcessLib::RichardsMechanics::MicroSolidVolumeFractionMode`. Optional;
default `"reference"`. A `<medium>` override defaults to the enclosing
block's mode when omitted.

Allowed values:
- `"reference"` (enum `MicroSolidVolumeFractionMode::Reference`): the
  active \f$n_S\f$ is held at the constant
  \ref ogs_file_param__prj__processes__process__RICHARDS_MECHANICS__potential_exchange__micro_solid_volume_fraction_reference "&lt;micro_solid_volume_fraction_reference&gt;"
  (clamped to at least \f$10^{-16}\f$) throughout the run
  (`computeActiveMicroSolidVolumeFraction()`,
  `computePreviousMicroSolidVolumeFraction()`,
  `ProcessLib/RichardsMechanics/RichardsMechanicsFEM-impl.h`). The
  associated tangent \f$dn_S/dn_l\f$ is `0`.
- `"current_porosity_split"` (enum
  `MicroSolidVolumeFractionMode::CurrentPorositySplit`): the active
  \f$n_S\f$ is instead the current aggregate solid fraction
  \f$1-n_l\f$ (clamped), re-evaluated from the current or previous
  micro-scale water content \f$n_l\f$ each time it is needed; the
  associated tangent \f$dn_S/dn_l=-1\f$ is folded into the micro-EOS and
  van-der-Waals-potential derivatives wherever this mode is selected.

When `"reference"` is selected, the full micro-disjoining-pressure
swelling-stress increment
(`computeReferenceMicroPorositySwellingStressIncrement()`) additionally
`OGS_FATAL`s unless `micro_solid_volume_fraction_reference > 0`.
