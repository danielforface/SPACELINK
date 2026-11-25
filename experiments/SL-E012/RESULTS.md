# E012 — energy-dependent mixture and scalar diffraction
2026-11-25. SIMULATED/DERIVED. Both Python 3.11 and 3.12 jobs completed successfully. CI https://github.com/danielforface/SPACELINK/actions/runs/45618327889. Revision cd22391f43df0e924926aff1885a143efd6c225e, exact JSON output. Pinned BS05 profile and Be7 spectrum retrieved/hash validated through E006 helpers. 92 energy samples and 183 midpoint-refined interpolation samples. Max Pee reduction disagreement .0001084784 (below .003), coverage 4.896856314e-9. Conditional p1 range .603217761–.604056168, p2 .373725948–.374562032, p3 .0222178834–.0222202070. 
| tungsten thickness | ideal recoil threshold MeV | conditional fractional count difference | full flux fractional count difference |
|---|---:|---:|---:|
| 1m | 0 | -2.72562e-13 | -1.32729e-21 |
| 1km | 0 | -2.01108e-7 | -9.79333e-16 |
| 1km | .25 | -2.00737e-7 | -9.77533e-16 |
| 1km | .5 | -1.98293e-7 | -9.65697e-16 |
At 1km, zero threshold, frozen p comparison -9.7978855e-16, energy dependence changes contrast magnitude ~0.0465%. Refined value -9.793329916e-16. Zero thickness residual ~5e-42 full relative due eigensolver roundoff, mathematical null exactly zero; not signal. No specialized external neutrino solver, measured detector response, residual coherence or actual solar spatial-energy correlation; threshold-only response. 
## Scalar grating model
Hypothetical Fourier component equal to illustrative tungsten electron potential, not a measured crystal. kappa=.000832103/m, Bragg angle .00359584 rad.
| normal length m | exact-on-resonance scalar model | angular first zero rad (Born) | angular averaged Born bound, ±.001rad |
|---:|---:|---:|---:|
| .001 | 6.92395e-13 | 7.19171e-10 | 2.48975e-19 |
| 1 | 6.92395e-7 | 7.19171e-13 | 2.48975e-16 |
| 10 | 6.92363e-5 | 7.19171e-14 | 2.48975e-15 |
The result is a scalar coupled-wave analogy with weak Born angular averaging, not a physical neutrino reflectivity or universal upper bound. It ignores spinor coupling, common neutral current, crystal structure, energy/angular joint distribution, crystal dimensions/thermal motion and receiver acceptance. Shows resonance versus acceptance tradeoff in toy model only. No claim actual steering, novelty or physical result. Next derive physical weak scattering amplitudes and compare S026, using structure factors, finite source and actual geometry.
