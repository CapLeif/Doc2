# Day 40 FL combined figure: FL01 approved reference

## Scope and evidence

This specification extracts the plot contract from Codex conversation `01a043a9-660a-7de0-90fb-d867e69dd4b8` and the local run record `KuroshioExt/FL_lambdaJPDF_20120101T060000/batch_NewFL/task_parameters.txt`. In that conversation, FL01 was rendered first, the user approved it, and the remaining figures were then generated. The batch record says 49 figures were written for IDs 2–50; it did not overwrite the approved FL01.

This is the Python Kuroshio LLC4320 FL01/batch figure family in this project. It is distinct from the MATLAB `εγθ_NewFL` contract already documented in the repository's `front.md` (which records a different data workflow and 600 dpi export). Do not merge their thresholds, sampling, or export resolution into a single universal rule.

The canonical executable implementation is included at `skills/doc2-day40-fl-figure/scripts/kuroshio_fl_theta_gamma_epsilon.py`, especially `plot_preview`, `draw_mode`, `arrow_samples`, `global_arrow_scale`, and `run_batch`. The figure is calculated from the original fields and extracted FL coordinates; it is not recreated from PNG pixels.

## Three-panel layout

- **Left, full height:** M² background in `summer`, cropped to a 14 km margin around the FL, a black smoothed centerline, white-centered sample markers, and four eigen-direction arrow groups. The panel uses equal x/y scale.
- **Upper right:** θ, signed γ, signed ε, two translucent `−λ/f` bands, and frontogenetic rate divided by f.
- **Lower right:** along-FL algebraic decomposition of the strain identity.
- Figure size is 15.0 × 9.2 inches, white background, exported at 220 dpi with tight bounds.

## Left-panel arrows and centerline

At each FL, nine sample positions are placed at equal arc-length fractions `0, 0.125, …, 1`. The plotted centerline is resampled every 0.25 km and smoothed in x(s), y(s) with a 2 km Gaussian scale; both endpoints are restored after smoothing. Smoothing is for display only. All field diagnostics remain sampled from the original grid.

Colors and directions:

| Mode | Color | Eigen-system |
|---|---|---|
| S− | light gray | smaller signed eigenvalue of symmetric strain tensor S |
| S+ | medium gray | larger signed eigenvalue of S |
| V− | blue | smaller real eigenvalue of Aᵀ |
| V+ | red | larger real eigenvalue of Aᵀ |

The display draws bidirectional axes; sign controls whether the pair is drawn inward from both ends (negative eigenvalue) or outward from the center (positive eigenvalue). Arrow half-length follows the direct square-root rule:

```text
half_length_km = 10 * sqrt(abs(lambda) / global_max_abs_lambda)
```

One global maximum is computed over the displayed FL set so lengths are comparable across figures. Do not replace this approved mapping with another scaling during reproduction.

## Upper-right diagnostics

The three rate values are:

```text
q = [-lambda_AT_plus/f, -lambda_AT_minus/f, frontogenetic_rate/f]
q_display = 90 + 30 * q / max(abs(q))
```

The two λ diagnostics are translucent bars spanning each point's exact arc-length bin edges, with baseline 90°, clipped to the 60–120° band. The V+ series is pale red, V− pale blue; white bin borders separate adjacent bars and alpha blending reveals overlap. The frontogenetic rate is a black/gray dashed curve mapped to the same scale, not a filled band. θ is red, signed γ green, and signed ε blue. The left axis spans −90° to 180°; a right secondary axis labels the diagnostic scale at 60°, 90°, and 120°.

## Lower-right strain identity

At each FL point, calculate:

```text
sigma = hypot(ux - vy, uy + vx)
zeta = vx - uy
delta = ux + vy
detLambda = ux*vy - uy*vx
det_part = -4*detLambda/sigma^2
zeta_part = zeta^2/sigma^2
div_part = delta^2/sigma^2
residual = det_part + zeta_part + div_part - 1
```

Require the largest absolute finite residual to be at most `1e-10` before saving the figure. Draw det contribution from zero (negative values below zero); stack ζ and δ terms algebraically above the positive det contribution. Expand y limits to include negative det and any positive sum above one. This makes the plotted contributions close according to the identity even where the determinant term is negative.

## Batch settings recorded for this run

```text
timestamp = 20120101T060000
top M2 percentile fraction = 10%
minimum FA major axis = 10 km
statistics FL length gate = strictly >20 km
raw candidates = 1361
raw FLs passing length gate = 398
continuous drawable segments = 397
plot_count = 50
rendered IDs after FL01 approval = FL02–FL50 (49 new files)
JPDF bin width = 0.5 degrees
arrow global maximum = 9.66244966279e-05 s^-1
arrow half-length at maximum = 10 km
```

These inventory values belong to the recorded input and run. Recalculate them for any changed fields, grid, selection, or parameters.

## Reproduction notes

The source script reads the Kuroshio LLC4320 binary fields at 960 × 931 × 40 and uses the first vertical level. It needs NumPy, SciPy, and Matplotlib. Example from the Skill folder:

```powershell
python scripts/kuroshio_fl_theta_gamma_epsilon.py `
  --root "D:\data\KuroshioExt" `
  --out "D:\results\day40-batch_NewFL" `
  --batch --plot-count 50 --plot-ids 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 `
  --top-percent 10 --min-major-axis-km 10 --stats-min-length-km 20 `
  --jpdf-bin-deg 0.5 --arrow-half-length-km 10
```

The code does not promise the same FL IDs or counts for a different field snapshot. Confirm the source timestamp and review one requested FL before a new large batch unless the user explicitly asks to render the full batch directly.

