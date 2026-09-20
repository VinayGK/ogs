`double`, optional, defaults to `1e-8` (a `<medium>` override that omits it
inherits the enclosing `<potential_exchange>` block's value). Must be
strictly greater than `0`, or the run stops with an
<tt>OGS_FATAL("local_jacobian_perturbation must be &gt; 0, ...")</tt>.

The source comment on the corresponding
`PotentialExchangeParameters` member describes it as the finite-difference
step for the implicit n_l(p_L) chain-rule derivative used in the
ScalarReferenceMassStorage local nonlinear-solve mode.

\attention As currently wired, the parsed value is stored on
`PotentialExchangeParameters` and range-checked at parse time, but is not
read back anywhere in the assembly code searched for this page. A dated
source comment at the ScalarReferenceStorage/ScalarReferenceMassStorage
tangent site records that the finite-difference dn_l/dp_L path this
parameter was written for ("perturbing the full coupled solve by
h ~ 1e-8*|p_L|") was replaced by an analytic tangent because it produced
catastrophic-cancellation noise at a dry initial condition; the tag survives
in the test fixtures
(`Tests/ProcessLib/RichardsMechanics/DSMMicroMacroSingleIntegrationPoint.cpp`)
but this page cannot point to a live production consumer of its value.
