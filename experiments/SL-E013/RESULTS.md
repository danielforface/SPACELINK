# E013 — ambient RF and parametric thermal budget
2026-11-25. DERIVED illustrative. Both Python 3.11 and 3.12 jobs completed successfully. CI https://github.com/danielforface/SPACELINK/actions/runs/45618327889. Revision cd22391f43df0e924926aff1885a143efd6c225e, 72 RF scenarios, 36 paired thermal scenarios x3 efficiencies, 4 asymmetric comparisons. JS independent central calculations agree to rounding.
## Radio
70 SFU, flat 100MHz at 2.8GHz, 1m² effective collect/transmit, 5m ideal receiver, D=225e9m, Tsys=100K: intercepted 7e-13W; optimistic reradiated receive 2.36830e-32W; AWGN reference 2.47474e-11 bit/s. Active 1W with same ideal aperture: 35.35336 bit/s. At 100m² collecting/transmitting area ambient 2.47474e-7 bit/s, active 3535.293 bit/s. Conditional lossless aperture model is favorable, includes all polarization, treats power as usable signal, no direct solar noise or extended incoherent source etendue; not a universal passive RF bound or physical observation. F10.7 scenarios, not actual day-by-day measured spectrum. No practical radio natural link established.
## Equal optical output 4W
320K radiator, emissivity .85, heat fraction .5 of unused solar power, collector 5kg/m², radiator 10kg/m²:
| assumed efficiency | collector m² | heat W | radiator m² | area-only mass kg |
|---|---:|---:|---:|---:|
| .0206 direct | .142671 | 95.0874 | .188145 | 2.594802 |
| .0396 optimistic direct | .074218 | 48.5051 | .095975 | 1.330834 |
| .09 PV+laser | .032656 | 20.2222 | .040013 | .563406 |
PV extra non-area mass break-even +2.0314 kg versus .0206 direct under these assumptions; +.76743kg vs .0396. Not actual component masses or beam quality. If direct .0206 rejects 99% of unused incident power without heating radiator and PV .09 absorbs all unused power, area-only mass .75098 vs .96353 kg: conditional reversal. Heat rejection fraction is hypothetical, not measured and not a claim real direct laser can achieve it. Need optical losses, heat flow, beam quality, component masses and temperatures before selection. No actual system mass advantage established.
