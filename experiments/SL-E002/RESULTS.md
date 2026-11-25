# SL-E002 — first screening run
Evidence type SIMULATED/DERIVED; S001 ideal channel arithmetic REPRODUCED, not physical BER reproduction. Code revision a8e50076c8deccdac1a2ccf5c2e2c67cd918781f; CI https://github.com/danielforface/SPACELINK/actions/runs/45379567953 ; both Python 3.11 and 3.12 jobs completed successfully, 10 tests OK (3.11 log https://github.com/danielforface/SPACELINK/actions/runs/45379567953/job/138394756382). Commands per protocol; standard library, deterministic seed 20261124, no external raw datasets. Parameters as preregistered.
## Results transcribed from 3.11 CI log
Independent dimensional conversion: 7.632466217868407e-14 eV per (g/cm3 * Ye), ~0.427% above rounded 7.6e-14. Keep rounded value for reproducibility; use precise source-audited constants in later model.
S001 ideal zero-background Z-channel: optimal prior 0.4050332072805001; capacity 0.36745428866347163 bits/pulse, grid 0.3674542872569458. Equal-prior uncoded BER 0.22242903311147055. Reported coded rate arithmetic 0.09881422924901184 bit/s (not independent BER measurement).
| Tungsten length m | Max leave probability | Persistent delta Pe | Two-flavor analytic | Two-flavor numerical |
|---:|---:|---:|---:|---:|
| 1 | 1.9552926541411082e-12 | -3.237241651739257e-13 | 1.9046881832513546e-12 | 1.9046881832513538e-12 |
| 100 | 1.9493213691988896e-8 | -3.2097011833513145e-9 | 1.904611839440082e-8 | 1.9046118394400816e-8 |
| 1000 | 1.817021878531626e-6 | -2.614530300383317e-7 | 1.8970651562011047e-6 | 1.897065156201105e-6 |
| 100000 | 0.0001513374108786314 | -2.1709722751295666e-5 | 0.00015824517915973298 | 0.00015824517915973663 |
100 seeded five-layer profiles, 300 input-state bound checks, max unitarity defect 4.000133557724439e-13, max positive bound excess 0.0. Analytic two-flavor agreement within preregistered tolerances; null vacuum and maximally mixed state checks pass. 
## Interpretation and limits
Standard matter can in this model change mass populations, so full post-material mass dephasing does not erase every possible change. However the persistent electron probability shift is tiny for these illustrative conditions. Fixed column bound limits ANY layering in this Hamiltonian regardless of vacuum gaps; no proof of global optimum, no experimental detection or sender-controlled natural flux link. Assumed mass probabilities, monoenergetic Be7, ideal complete dephasing, no detector response, no finite source geometry. Do not combine Pe shift directly with pp 134 events/day; Be7 fitted rate ~48.3/100t/day, flavor-dependent cross section required. External solver and independent human review outstanding. No claim of novelty.
