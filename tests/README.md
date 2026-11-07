# tests — unit tests with a known answer

Run from the repository root:
```matlab
runtests('tests')
```
21 tests, all passing on MATLAB R2025a (2026-11-07). Each file adds `../src` to the path.

| file | test | checks |
|---|---|---|
| `test_core.m` | `test_partitionOfUnity` | Bernstein polynomials sum to 1 for n = 0..8 |
| | `test_classFunctionEllipse` | $N_1 = N_2 = 0.5$, $S = 1$ is a circle |
| | `test_leadingEdgeRadius` | $A_0 = \sqrt{2 r_{LE}/c}$ |
| | `test_fitCstCurveExact` | a CST curve is fitted exactly |
| | `test_spanwiseBasisContinuity` | piecewise basis continuous at the breaks, end interpolation, partition of unity |
| | `test_fitSpanwiseTooFewStations` | clear error when a segment has too few stations |
| `test_airfoil.m` | `test_roundTrip` | CST → airfoil → CST, coefficients to 1e-10 |
| | `test_normalizeRotatedScaled` | chord, twist, LE recovered from a scaled, rotated, moved airfoil |
| | `test_fileAndOrder` | Selig file, reversed point order, thickness of NACA 0012 |
| | `test_fitAccuracy` | NACA 0012 at n = 6 within 2e-4, NACA 2412 at n = 8 |
| | `test_lednicer` | Lednicer file gives the same airfoil |
| `test_wing.m` | `test_roundTrip` | cranked wing CST → wing → CST |
| | `test_planformGeometry` | area, aspect ratio, sweep of a tapered wing |
| | `test_fromSections` | planform from 3D sections |
| | `test_breakNeeded` | the kink is exact with the break, > 0.1 m chord error without |
| | `test_normalizeSectionsKeepsGeometry` | normalizing the section airfoils does not move the 3D sections |
| `test_nacelle_fuselage.m` | `test_nacelleRoundTrip` | nacelle CST → nacelle → CST |
| | `test_nacelleSmoothAtSymmetryPlane` | zero derivative in $\varphi$ at crown and keel |
| | `test_nacelleAxisymmetric` | one profile → constant radius around |
| | `test_fuselageRoundTrip` | fuselage CST → fuselage → CST |
| | `test_fuselageFromSections` | lines and lobe exponents recovered from cross sections |
