# MITgcm LLC4320 processing flow in this project

## Scope

This guide records the actual Python processing paths found in the project. “LLC4320” is not one universal file contract here: the Kuroshio line workflow reads a fixed 960 × 931 × 40 big-endian binary snapshot, West Atlantic uses NetCDF, and the generic regional surface script expects 2-D fields. Keep these paths separate and inspect each input before running.

## Main processing paths

| Goal | Source | Input contract | Result |
|---|---|---|---|
| Surface θ, λ−, Ro, σ/|f|, divergence/f maps and statistics | `skills/mitgcm-llc4320-processing/scripts/process_surface_diagnostics.py` | One region folder with 2-D `XC`, `YC`, `U`, `V`, `hFacC`, `DYC`, optionally `DXC`, `AngleCS`, `AngleSN`; big-endian float32, Fortran order | Compressed NPZ, field PNGs, diagnostic suite, statistics PNG/TXT |
| Full Kuroshio FA/FL and along-front figure/statistics/JPDF | `skills/mitgcm-llc4320-processing/scripts/kuroshio_fl_theta_gamma_epsilon.py` | Fixed LLC Kuroshio names and shapes at 960 × 931 × 40; big-endian float32, Fortran order; timestamp is `20120101T060000`; takes level index 0 | FL preview or batch PNGs, point/summary tables, angle JPDF, task parameters |
| Regional FA/FL angle JPDF | `skills/mitgcm-llc4320-processing/scripts/regional_fl_jpdf.py` | Binary mode uses 2-D region files; West Atlantic mode reads the configured NetCDF file | Inventory, point CSV, angle JPDF CSV/PNGs, task parameters |
| Three-region velocity-gradient λ-JPDF | `skills/mitgcm-llc4320-processing/scripts/plot_three_region_uvgrid_lambda_jpdf.py` | Precomputed surface NPZ for KE/SCS; West Atlantic NetCDF for local-f recalculation | Three 600 dpi λ-JPDF PNGs and one summary CSV |
| Eigenmode classification fields and conditional λ distributions | `skills/mitgcm-llc4320-processing/scripts/kuroshio_mode_classification.py` | Inputs loaded by the Kuroshio and regional workflow modules; inspect CLI/variables before running | Classification NPZ, figures, tables |

## Kuroshio LLC4320 full-field FL path

### 1. Read the snapshot

The script at `skills/mitgcm-llc4320-processing/scripts/kuroshio_fl_theta_gamma_epsilon.py` fixes `NX=960`, `NY=931`, `NZ=40`, `STAMP=20120101T060000`. It reads XG/YG, DXC/DYC, hFacC and Theta/Salt/U/V files with shape `(960,931)` or `(960,931,40)` as named, `dtype='>f4'`, and `reshape(order='F')`. Theta, Salt, U, and V are sliced at `[:,:,0]`. Data are not dynamically detected by timestamps; change the code/configuration only after validating a new snapshot's file dimensions and variable layout.

Wet/finite screening requires `hFacC>0`, legal finite coordinates, finite scalar and velocity values, and `|U|,|V|<20`. U is collocated by averaging adjacent values along axis 0; V along axis 1. The buoyancy definition is `b=g(αT·Theta−βS·Salt)`, with `g=9.81`, `αT=2e−4`, and `βS=7.6e−4`. Spatial derivatives are `np.gradient(field, axis)/DXC` or `/DYC`. The horizontal axes are projected to km using a local equirectangular approximation at median latitude. The script computes `f=2Ωsin(lat)`, with `Ω=7.2921159e−5 s−1`, and `M²=hypot(bx,by)`.

### 2. Find high-gradient regions and front lines

The default top-10% M² mask uses the finite-field 90th percentile and 8-neighbor labeling. `M²` is nearest-filled and smoothed with a one-grid Gaussian for ridge derivatives. The ridge correction uses the buoyancy-gradient normal, negative normal curvature, and a sub-grid offset no larger than half a grid diagonal. Candidate ridge pixels form a tangent-aligned 8-neighbor graph (default maximum departure 45°); branch nodes are removed before connected paths are summarized. Candidates are ordered by length and then filtered by requested thresholds. These are code-defined operational parameters; recalculate thresholds on the actual finite input domain.

### 3. Compute along-line quantities

The script interpolates grid fields onto FL pixels, constructs `A=[[ux,uy],[vx,vy]]`, its symmetric part `S=(A+Aᵀ)/2`, and the eigenvectors/eigenvalues needed for θ, α, β, signed γ/ε, normalized λ values and frontogenetic rate. It also evaluates `ζ²/σ² + δ²/σ² − 4det(A)/σ² = 1`; the figure path stops if the largest finite residual exceeds `1e−10`.

### 4. Render and summarize

`--batch` writes combined figures and tables under the selected output root. The code uses the requested top percentile, 10 km minimum FA major axis, `>20 km` default statistics length, 0.5° default JPDF bins, 50-figure limit, 10 km maximum arrow half-length, and a single maximum-eigenvalue scale for the displayed atlas. The entry point can accept `--root`, `--out`, `--top-percent`, `--min-major-axis-km`, `--stats-min-length-km`, `--plot-count`, `--plot-ids`, and `--jpdf-bin-deg`. Check the actual run's `task_parameters.txt`; never report another snapshot's counts as current.

Example:

```powershell
python skills/mitgcm-llc4320-processing/scripts/kuroshio_fl_theta_gamma_epsilon.py `
  --root "D:\data\KuroshioExt" --out "D:\results\KuroshioExt-Day40" `
  --batch --plot-count 50 --top-percent 10 --min-major-axis-km 10 `
  --stats-min-length-km 20 --jpdf-bin-deg 0.5 --arrow-half-length-km 10
```

## Other LLC4320/region routes

### Generic 2-D surface diagnostics

`skills/mitgcm-llc4320-processing/scripts/process_surface_diagnostics.py` reads 2-D binary fields and collocates U and V to tracer points before computing `ζ=vx−uy`, `σ=hypot(ux−vy,vx+uy)`, `δ=ux+vy`, `λ−=0.5(δ−sqrt(σ²−ζ²))`, and `θ=acos(−ζ/σ)`. It masks non-real eigenvalue points (`σ²−ζ²≤0`) and small `|f|≤1e−8`. If DXC is absent, it computes a local curvilinear geographic gradient from longitude/latitude. If AngleCS/AngleSN are both present, it rotates model-grid velocities. This script is not a direct replacement for the fixed 3-D Kuroshio FL loader.

### Regional FA/FL JPDF

`skills/mitgcm-llc4320-processing/scripts/regional_fl_jpdf.py` provides a `binary` mode for supported 2-D region files and a `westatlantic-nc` mode for the named LLC4320 NetCDF layout. In the latter, it trims the extra boundary row/column, collocates face velocities, rotates with AngleCS/AngleSN, computes buoyancy and geographic gradients, and uses local latitude for f. If the binary scalar path has no usable salinity and the explicit temperature-gradient option is chosen, the scalar is a temperature-gradient proxy, not complete buoyancy; preserve that distinction in output names and interpretation.

### Three-region λ-JPDF

`skills/mitgcm-llc4320-processing/scripts/plot_three_region_uvgrid_lambda_jpdf.py` applies KE bounds 155–165°E and 36–40°N, South China Sea bounds 115–125°E and 17–23°N, and full West Atlantic. The first two consume prepared NPZ diagnostics; the West Atlantic branch recomputes velocity-gradient fields and f from its NetCDF. It maps the real-eigenvalue region to negative `−sqrt((σ/|f|)²−Ro²)` and the complex-pair region to positive `sqrt(Ro²−(σ/|f|)²)`, then estimates a 0.05-bin PDF over `[-10,10]²`. Set `MITGCM_PROJECT_ROOT` to the local project/data root before running this script; it does not infer a data directory from an installed Skill path.

## Validation and portability

- Record input paths, timestamp, file dimensions, byte order, reshape order, level index, grid staggering, masks, and parameters for every run.
- On a new grid/snapshot, compare at least one derived field with a trusted existing diagnostic before a multi-region batch.
- Separate the fixed 3-D Kuroshio FL path from generic 2-D surface diagnostics and West Atlantic NetCDF loading; they are not interchangeable readers.
- Use local data paths on each computer. This repository intentionally contains scripts and instructions only, not raw MITgcm datasets or generated scientific results.
- Python packages required by the included scripts are listed in `skills/mitgcm-llc4320-processing/scripts/requirements.txt`.

