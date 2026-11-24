# SL-E001 — bootstrap screening run
Evidence: DERIVED screening model; execution verified in CI. Not an observed signal, not a detector prediction, not scientific certification. Source constants and bound still pending independent review.
Code commit: 69bed0219ea95235cd2dac219976c0a410acb341
CI: https://github.com/danielforface/SPACELINK/actions/runs/45377180231
Python 3.11 job: https://github.com/danielforface/SPACELINK/actions/runs/45377180231/job/138386454203
Python 3.12 job: https://github.com/danielforface/SPACELINK/actions/runs/45377180231/job/138386454028
Both jobs successful; test command and calculation command succeeded. Python 3.11 log records 5 tests OK. No external dataset used. Local git clone failed DNS resolution; no local test run claimed.
Commands: python -m unittest discover -s tests -v; python experiments/SL-E001/bounds.py.
Values transcribed from Python 3.11 CI log:
| Material | Length m | Phase rad | Transition upper bound |
|---|---:|---:|---:|
| water | 1 | 2.1375688167171692e-7 | 1.1423001115504054e-14 |
| water | 100 | 2.137568816717169e-5 | 1.1423001115069144e-10 |
| water | 1000 | 0.00021375688167171692 | 1.1423001072009114e-8 |
| water | 100000 | 0.02137568816717169 | 0.00011422566172280335 |
| tungsten | 1 | 2.9921002438817353e-6 | 2.2381659673576154e-12 |
| tungsten | 100 | 0.0002992100243881735 | 2.2381659506613283e-8 |
| tungsten | 1000 | 0.0029921002438817353 | 2.2381642975641506e-6 |
| tungsten | 100000 | 0.29921002438817357 | 0.02221517761847063 |
Ideal toy integration days at 1% fractional detector modulation: 1865.6716417910447. This assumes known baseline, no backgrounds, approximate 134 events/day and 100 t at 1 AU; unrelated to actual modulation achievable by material. The 5-sigma screening approximation is not a communication rate or BER.
No source audit, independent solver or uncertainty propagation completed yet. Follow issue #1 for scientific verification and protocol SL-E002.
