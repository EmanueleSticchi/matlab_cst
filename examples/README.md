# examples — build and fit every shape

**Created:** 2026-11-07 · **Status:** ready · **MATLAB:** R2025a (the `build_*` live scripts need R2025a or later)

Run every script from this folder: they add `../src` to the path with a relative path.

## Build: live scripts with sliders
`build_*.m` are **live scripts** saved in the plain-text Live Code format (R2025a+): open them in the Live Editor (double click). Every CST parameter is a slider; moving it re-runs the second section, which builds the shape with the `CST2*` function and plots it. The default values come from fits of real shapes, written in the text at the top of each script.

| script | CST parameters on the sliders | plots |
|---|---|---|
| `build_airfoil.m` | `Au`, `Al` (order 4), `dzu`, `dzl`, `N1`, `N2` | airfoil with its component airfoils; shape functions $S_u$, $-S_l$ and their Bernstein terms |
| `build_wing.m` | span, break position, chord / Xle / Zle / twist at root, break and tip, root and tip `Au`, `Al` (order 3, `Al(1) = -Au(1)`) | 3D wing (both halves), planform, root / break / tip airfoils |
| `build_nacelle.m` | `Rle`, `c`, `Xle`, twist and `Au`, `Al` on the crown and keel lines | 3D nacelle, crown / side / keel profiles |
| `build_fuselage.m` | `L`, `N1`, `N2`, `Ahu`, `Ahl`, `Aw` (order 4), tail values, `Bzm` (order 2), `NU`, `NL` | 3D fuselage, side view of the lines, largest section |

In releases before R2025a the files still run as plain scripts (the slider lines are ordinary assignments).

## Fit: from the physical domain to the CST
`fit_*.m` follow the script skeleton (DATA / PROCESS / PLOT / FUNCTIONS): everything to change is in DATA.

| script | input geometry | what it shows |
|---|---|---|
| `fit_airfoil.m` | `data/naca0012.dat` (any Selig or Lednicer file) | fit, residuals on the two sides, convergence with the Bernstein order (Kulfan Fig. 9–11) |
| `fit_wing.m` | cranked wing, 7 NACA sections, kink at $\eta = 0.35$ (`Wing3D_class`) | fit with and without the break: planform error 1e-5 m vs 0.35 m |
| `fit_nacelle.m` | 5 cowl profiles crown → keel, scarfed and drooped (`Nacelle_class`); an axisymmetric nacelle | input vs CST profiles, 3D |
| `fit_fuselage.m` | 60 synthetic cross sections (point clouds), nose / cylinder / upswept tail (`Fuselage_class.fromSections`) | lines input vs CST, convergence with the order |

With `saveflag = true`, `fit_airfoil.m` and `fit_wing.m` write `*_cst.mat` (the CST struct) and `*_cst.dat` (Selig airfoil / Tecplot wing) here; they are git-ignored.

## Data
| file | source |
|---|---|
| `data/naca0012.dat` | NACA 0012, 120 panels, half cosine spacing (copied from WiMOT `examples/naca0012.dat`) |

Other NACA 4-digit airfoils are generated with `Airfoil_class.naca4("2412")`.

## Results (2026-11-07, MATLAB R2025a)
| script | key numbers |
|---|---|
| `fit_airfoil.m` | NACA 0012, n = 6: z/c error rms 3.8e-5, max 1.4e-4 |
| `fit_wing.m` | with break: chord error 2e-5 m, Xle 1e-5 m; without: 0.36 m, 0.19 m; airfoil rms ≤ 8e-4 (NACA catalogue steps along the span) |
| `fit_nacelle.m` | airfoil rms ≤ 7e-4, Rle error 5e-4 m |
| `fit_fuselage.m` | max error on hU, hL, w: 3.5 cm (n = 4) → 1.2 cm (n = 12), 2 m radius |
