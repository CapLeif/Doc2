# Doc2 Core Conventions

Compact reference for the shared variables, partitions, algorithms, and figure conventions used in the Doc2 flow and frontal analyses.

## Project-wide conventions

- The horizontal velocity-gradient matrix is \(\Lambda=\nabla_h\mathbf{u}=\begin{bmatrix}u_x&u_y\\v_x&v_y\end{bmatrix}\). Its transpose is written \(\Lambda^\mathsf{T}\). Do not use \(\Lambda\) and \(\Lambda^\mathsf{T}\) interchangeably: their eigenvectors differ.
- The strain-rate tensor is \(S=(\Lambda+\Lambda^\mathsf{T})/2\). A negative signed eigenvalue denotes compression; a positive signed eigenvalue denotes extension.
- \(\alpha=\beta\) geometrically in the earlier paired naming; this project writes that angle only as \(\alpha\), with no separate \(\beta\) variable.
- Unless stated otherwise, angles are in degrees, spatial coordinates and along-front distances are in km, velocity-gradient quantities are in s\(^{-1}\), and Coriolis-normalized quantities are dimensionless.
- Use \(|f|\) where a nonnegative strain magnitude is normalized; retain signed \(f\) where the diagnostic defines a signed ratio such as \(\zeta/f\), \(\Delta/f\), or \(\lambda_-/f\). Each axis label must show which convention applies.
- Missing, nonfinite, zero-denominator, complex-eigenvalue, and degenerate-eigenvector cases are excluded only as specified by the individual diagnostic. Do not silently replace them with zero.

## Abbreviations

| Short name | Meaning |
|---|---|
| VGL | Compact index for the project's named flow-partition labels in [part.md](part.md). No expansion of the letters is asserted. |
| JPDF | Joint probability density function |
| FA | Frontal area |
| FL | Frontal line |
| SD | Strain dominated |
| VD | Vorticity dominated |
| AVD / CVD | Anticyclonic / cyclonic vorticity dominated |
| M2 | Project input field used to identify high-gradient frontal areas; its upstream derivation is not redefined by the plotting/extraction workflow. |
| Ro | Project input Rossby-number field; use the saved field for masks unless a diagnostic explicitly recomputes it. |

## Documents

- [vars.md](vars.md) — variables, equations, units, and angle definitions.
- [part.md](part.md) — VGL partition names, criteria, and numeric identifiers.
- [jpdf.md](jpdf.md) — lambda-JPDF coordinates, histogram, overlays, and labels.
- [field.md](field.md) — parameter-field masks, positive-only ranges, and map style.
- [front.md](front.md) — frontal-area/ridge extraction and FL eigenvector-arrow figures.

## Common source basis and applicability

These conventions were extracted from the project manuscript and MATLAB sources, principally `Field_on_1234.m`, `M6_sigma_delta_zeta_day40/src/compute_frontogenesis_diagnostics.m`, `M6_sigma_delta_zeta_day40/src/front_compute_fl_angles.m`, `M6_sigma_delta_zeta_day40/src/front_extract_fullfield_fas.m`, `M6_sigma_delta_zeta_day40/src/front_extract_ridges.m`, and `M6_sigma_delta_zeta_day40/scripts/run_day40_newfl_theta_gamma_epsilon.m`. The manuscript partition table and Figure 2 description supply the VGL regime names and inequalities.

The sources include figure- and dataset-specific configurations. Their exact values are identified as examples or local settings in the relevant document; they are not universal physical constants. The FL01 reference is `FL_lambdaJPDF/εγθ_NewFL/combined_figures/FL01_arrows_theta_gamma_epsilon.png`; its drawing contract is read from the generating MATLAB code.

## Release scope

This repository contains definitions and conventions only. It excludes source data, intermediate data, generated figure collections, and the full analysis code. Each topic file contains only topic-specific rules; shared project context and source provenance belong here.

