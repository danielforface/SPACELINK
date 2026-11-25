# E008 — variable amplitude and resolution results
2026-11-25, DERIVED/SIMULATED. Code 6eb31cd5aeca3f32481cb7048924a6b04d712f1b, run https://github.com/danielforface/SPACELINK/actions/runs/45588052309, Python 3.12 job 139114648189 logs 40 tests OK, E008 and E009 passed. Raw outputs in results-ci.json. Initial code d2fe7e2 run 45587891891 passed 39 tests, unchanged E008 physics. Both corrected CI jobs completed successfully. No physical measurement.
## Variable-amplitude bound
Fixed BS05 aperture-conditioned initial mixture and tabulated Be7 piecewise-linear normalized spectrum, uniform tungsten. Derivation in PROTOCOL. 
| thickness m | amplitude norm bound | derivative bound /eV | sampled maximum off-diagonal | sampled maximum derivative /eV | 21 entry bound at .524 AU |
|---:|---:|---:|---:|---:|---:|
| 1 | 3.61348715e-6 | 3.09761138e-14 | 3.10298046e-7 | 9.10714386e-16 | 1.79597881e-10 |
| 1000 | .00361348715 | 3.09761138e-8 | .000309683164 | 2.18348013e-10 | 1.81165573e-7 |
31 bounds 5.38793644e-12 and 5.43496718e-9; 32 bounds 5.55457365e-12 and 5.60305895e-9. Finite-difference and short-baseline quadrature checks passed; analytical bound is conservative; sampled maxima do not establish global extrema. Fixed mixture, no energy dependence in initial populations or potential, no detector response, no real measurement, not universal no-go or capacity bound. 
## Phase resource scales at .524 AU (78389284246.8m)
| pair | full energy bin eV (phase span <=1 rad) | position variation m / rad | separation m | illustrative Gaussian sigma_x for half overlap m |
|---|---:|---:|---:|---:|
| 21 | .04987845987 | 4535.8895 | 3.95615624e-6 | 1.68002487e-6 |
| 31 | .00149635380 | 136.0767 | .000131871875 | .0000560008290 |
| 32 | .00154263278 | 140.2852 | .000127915719 | .0000543208042 |
Independent JS computation confirmed width/position/separation; Gaussian half-contrast energy sigma values in JSON. Protocol mentioned e^-1 overlap as well; output only includes half overlap (e^-1 follows analytically sigma_x=delta_x/sqrt(8)); no actual solar or modulator packet size established. Do not multiply independent wave-packet and energy damping factors without consistent detection model (S041). Still lack true detector energy-response and energy-dependent production mixture integration, wave-packet QFT and external physics validation. 
