# Doc2 Reference Center

Cloud reference for Doc2 variables, algorithms, MATLAB implementations, and figure conventions. Use this repository as the project-wide source of truth when planning or producing related diagnostics and figures.

## Project-wide conventions

- The horizontal velocity-gradient matrix is `Λ = ∇h u = [[u_x,u_y],[v_x,v_y]]`. Its transpose is `Λᵀ`; they are not interchangeable.
- The strain tensor is `S = (Λ + Λᵀ)/2`. Negative signed eigenvalues denote compression; positive signed eigenvalues denote extension.
- `α = β` geometrically in earlier paired naming. This project retains `α` only and does not define a separate `β` variable.
- Unless a topic says otherwise, angles are in degrees, distances in km, velocity-gradient quantities in s⁻¹, and Coriolis-normalized quantities are dimensionless.
- Use `|f|` when normalizing a nonnegative strain magnitude. Preserve signed `f` for signed ratios such as `ζ/f`, `Δ/f`, and `λ₋/f`; label axes accordingly.
- Missing, nonfinite, zero-denominator, complex-eigenvalue, and degenerate-eigenvector cases follow each diagnostic's stated validity rule; do not replace them silently with zero.
- The authoritative domain/partition mask is supplied by the calling workflow. Do not infer VGL meanings from numeric category IDs unless the mapping is confirmed in `part.md`.

## Abbreviations

| Short name | Meaning |
|---|---|
| VGL | Project's named flow-partition labels; no expansion is asserted here. |
| JPDF | Joint probability density function |
| FA / FL | Frontal area / frontal line |
| SD / VD | Strain dominated / vorticity dominated |
| AVD / CVD | Anticyclonic / cyclonic vorticity dominated |
| M2 | Input field used by the frontal-area workflow; upstream derivation is outside this plotting/extraction contract. |
| Ro | Rossby-number input field; use the saved field for masks unless a diagnostic explicitly recomputes it. |

## Reference documents

- [vars.md](vars.md): strict variables, equations, units, and relationships.
- [part.md](part.md): VGL names, criteria, and numeric identifiers; unverified mappings are marked.
- [jpdf.md](jpdf.md): λ-JPDF coordinates, binning, normalization, contours, and labels.
- [field.md](field.md): parameter-field masks, value ranges, colorbars, and map display.
- [front.md](front.md): frontal-area/ridge extraction and FL feature-direction arrows.

## Core MATLAB code

The `code/` tree contains selected reusable implementations. Function filenames retain their MATLAB entry-point names; the theme folders make the collection easier to browse. These modules accept project arrays/structs from the caller and do not include data, machine-specific paths, or dataset-specific run scripts.

| Theme | Main implementation | Local source |
|---|---|---|
| Variables/diagnostics | `code/vars/compute_frontogenesis_diagnostics.m`, `compute_theta_from_zeta_sigma.m`, `compute_m6_residual.m` | `M6_sigma_delta_zeta_day40/src/` |
| Partition/flow diagnostics | `code/part/front_compute_diagnostics.m`, `front_compute_at_lambda_minus.m`, `front_compute_fl_angles.m`, `front_compute_fl_lambda_jpdf.m`, `front_classify_structure.m` | `M6_sigma_delta_zeta_day40/src/` |
| λ-JPDF | `code/jpdf/compute_lambda_conditional_means.m`, `front_plot_fl_lambda_jpdf.m`, `plot_lambda_frontogenesis_conditionals.m` | `M6_sigma_delta_zeta_day40/src/` |
| Parameter fields | `code/field/field_filter_mask.m`, `plot_frontogenesis_fields.m`, `plot_m6_residual_field_combined.m` | `field_filter_mask.m` is the reusable form of mask logic in `Filaments/plot_day40_density_ro_lambda_filters.m`; other files from `M6_sigma_delta_zeta_day40/src/` |
| Front extraction/arrows | `code/front/front_extract_fullfield_fas.m`, `front_extract_ridges.m`, `front_sample_eigenvectors.m`, `front_plan_fullfield_eigen_samples.m`, `front_plot_fullfield_eigenvectors.m`, `front_plot_single_fl_eigenvectors.m`, `front_plot_curated_single_fl.m` | `M6_sigma_delta_zeta_day40/src/` |

`field_filter_mask` takes an authoritative base mask from the caller and applies the shared finite-value, `|Ro| ≥ 1`, and `λ₋/f < -1` filters. It deliberately does not construct a non-SD mask from numeric VGL IDs. Positive-only display limits and colorbar settings remain governed by `field.md` and should be supplied by the plotting workflow.

The selected functions are reference implementations, not a standalone application: use the documented input fields, MATLAB version/toolboxes, and optional plotting utilities required by each function. The full Day40 orchestration scripts remain local because they contain dataset-specific paths/configuration; their reusable conventions are recorded in the topic documents.

## Provenance and applicability

The reference is distilled from the Doc2 manuscript and local MATLAB sources, principally `Field_on_1234.m`, `Z_JPDF_ROMS.m`, `M6_sigma_delta_zeta_day40/src/compute_frontogenesis_diagnostics.m`, `front_compute_fl_angles.m`, `front_extract_fullfield_fas.m`, `front_extract_ridges.m`, and the local NewFL/field plotting workflows. Topic documents identify their source basis. Figure- and dataset-specific values are examples unless explicitly marked as project-wide rules.

The FL arrow contract uses the local generating code associated with `FL_lambdaJPDF/εγθ_NewFL/combined_figures/FL01_arrows_theta_gamma_epsilon.png`; the image itself is not redistributed in this repository.

## Repository contents and maintenance

`SKILL.md` routes future Doc2 figure and diagnostics tasks to the relevant documents and code. Keep shared symbols, units, provenance, and applicability here in this README; keep topic rules in their one topic document. When changing a convention, update its topic document, this shared index if needed, and the corresponding core implementation together.

This public repository excludes source/intermediate data, generated figure batches, machine-specific paths, and full dataset-specific run scripts. Version history is the change record; no separate release version is asserted yet.

