# E011 — full count versus binary receiver
2026-11-25. SIMULATED ideal independent Poisson gate. Both Python 3.11 and 3.12 jobs completed successfully. CI https://github.com/danielforface/SPACELINK/actions/runs/45618327889. Revision cd22391f43df0e924926aff1885a143efd6c225e; 46 tests OK in job 139219468035. 120 scenarios; max numerical retained mass error 2.5161e-10 (not a rigorous floating-point tail certificate). 100 seeded checks and recurrence comparisons passed. Initial CI 45618220059 failed in independent recurrence due subnormal denominator; corrected before rerun, protocol amendment and original CI retained.
## Same photons, same prior and gate
10nm passive mirror, p=.5, gate 1ms, guard 1us, 58.53360665 signal counts/s:
| Synthetic radiance ratio | Background counts/s | Full-count bits/s | Binary click bits/s | Uncoded count BER |
|---|---:|---:|---:|---:|
| 1e-10 | 162.7465 | 3.220089 | 2.929211 | .475843 |
| 1e-7 | 62846.4793 | .009816913 | 3.046993e-28 | .498528 |
| 1e-4 | 62746579.3 | 9.837184e-6 | numerical zero | .499953 |
Full count restores information discarded by click saturation; this follows data processing and is not novel. Middle scenario ratio ~3.22e25 is relative to an intentionally saturated baseline; absolute rate remains below .01 bit/s. BER near .5 per raw gate; asymptotic MI not a demonstrated coded link. At active 1nm, bg 1e-7, 100ns gate, 20ns guard, p=.5, both direct solar and assumed PV full count ~8.333333Mbps vs click ~8.301098Mbps; one-bit-per-gate ceiling. No hardware count throughput, dead time, afterpulsing or analog response; high background rates in table may be impossible for a real detector. No link margin or physical observation. Next finite detector count/analog dynamic range, coding and E010 bench evaluation after equipment/safety review.
