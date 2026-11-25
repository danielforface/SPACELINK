# SL-E003 — geometry/detection screen results
Evidence: SIMULATED/DERIVED under preregistered simplified assumptions. Not a physical observation or scientific certification. Code abd012474d1f9088e80cfce638d72a78c9f42319; CI https://github.com/danielforface/SPACELINK/actions/runs/45380593988 both Python 3.11/3.12 jobs passed; 15 tests OK in Python 3.11 log https://github.com/danielforface/SPACELINK/actions/runs/45380593988/job/138398295676. NumPy 1.26.4; exact numerical output transcribed from CI into results-ci.json. Commands per protocol.
## Independent numerical oracle
112 cases; maximum absolute probability disagreement 4.318767565791859e-13; maximum relative off-diagonal probability error 6.689782061641836e-11 (where reference P>1e-10). Within preregistered tolerances. Independent numerical algorithm (LAPACK eigh versus Taylor scaling/squaring), shared Hamiltonian/assumptions, not an independent researcher or specialized neutrino solver. No verification of solar initial population model.
## Receiver
Integrated tree-level electron-to-other flavor scattering cross-section ratio 4.76123463182996. Analytic integral agrees with independent trapezoid quadrature in tests. Be7 central fitted rate proxy 48.3/day/100t at 1 AU => 20.795 approx at 1.524 AU; no claim of raw event availability or realistic detection calibration.
## Link screen
100t receiver; exact alignment; uniform projected disk effective radius .1 solar radii; no backgrounds. State 0 vacuum, state 1 ordinary tungsten, complete post-material mass dephasing. Rates are upper bounds for binary total-count Poisson channel only, not predicted realizable throughput.
| Thickness m | Radius m | Covered fraction | Abs relative count change (simulation) | Rate upper bits/day |
|---:|---:|---:|---:|---:|
| 1 | 1 | 1.7476857828e-15 | 6.9227715685e-28 | 2.9554796382e-52 |
| 1000 | 1 | 1.7476857828e-15 | 5.5911167518e-22 | 2.9554752283e-40 |
| 1000 | 1000 | 1.7476857828e-9 | 5.5911167518e-16 | 2.9554752283e-28 |
| 1000 | 100000 | 1.7476857828e-5 | 5.5911167518e-12 | 2.9554752283e-20 |
| 1000 | 10000000 | .1747685783 | 5.5911167518e-8 | 2.9554752283e-12 |
| 1000 | full coverage hypothetical | 1 | 3.1991544514e-7 | 9.6761059915e-11 |
At 1km thickness/radius, 100 background events/day changes bound to 5.0880622563e-29; proxy-rate sensitivity 46.1–50.5 yields 2.8208573e-28 to 3.0900931e-28 without background. This is scenario sensitivity, not uncertainty interval.
## Interpretation
Severe negative result for these assumptions. Total-count dilution f gives f² information-rate penalty, for fixed other parameters. Perfect direction tagging could remove unmodulated foreground but sacrifices event count (ideal upper-bound scaling f, not f²); not implemented/verified as a practical receiver. No geometry-independent no-go theorem. No actual source distribution ingested: IAS file header inspected, runtime DNS blocked. Uniform disk is illustrative, not rigorous bound on actual solar geometry. Recoil-energy and angular response, finite detector, full solar spectrum and nonaligned geometry remain. See literature/SL-P002 and protocol. Next: test ideal tagged-event bound, ingest actual shell distribution if available, and prioritize other interaction mechanisms before expensive profile optimization. No novelty claim.
