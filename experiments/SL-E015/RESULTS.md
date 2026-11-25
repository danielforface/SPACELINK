# E015 — isolated nuclear scattering angular screen
2026-11-25. DERIVED/SIMULATED standard leading vector neutral-current approximation with illustrative uniform spherical weak charge. Both Python 3.11 and 3.12 jobs completed successfully. CI https://github.com/danielforface/SPACELINK/actions/runs/45624762394. Exact rows in JSON; 9 material-energy cases, 4 angular cones per case, 128/256 quadrature agreement, analytical point-nucleus cone checks. Source S060 Eq8–11. 
W-184 at .862MeV: point cross section 3.4259465454e-41 cm², uniform-sphere 3.4251345272e-41 cm², ratio .9997629799; 1km optical depth 2.1635554094e-13. At 10 and 50MeV form-factor integrated ratios .9687417445 and .4975531406. Not a measured W form factor. Heavy nucleus approximation, spin zero, no axial/electron corrections. 
| half cone rad | fraction of scattered events forward | backward | probability per incident particle for 1km W, forward |
|---:|---:|---:|---:|
| 1e-8 | 5.00119e-17 | 6.24704e-34 | 1.08203e-29 |
| 1e-4 | 5.00119e-9 | 6.24704e-18 | 1.08203e-21 |
| .01 | 5.00108e-5 | 6.24693e-10 | 1.08201e-17 |
| .1 | .00499077 | 6.23664e-6 | 1.07978e-15 |
Forward fraction scales alpha², backward alpha⁴ for small angles in leading helicity distribution. Backward cone is around reverse of incident direction; these are not arbitrary-angle steered beams. Full reversal suppressed by 1+cos theta. Nuclear coherence barely reduces total cross section at solar energy; no enhancement of macroscopic crystal diffraction demonstrated. Independent scatter geometry, not coherent forward propagation or lattice Bragg calculation. Specialized neutrino propagation oracle still not run, actual crystal weak structure factor and wave packet model outstanding. No physical result.
