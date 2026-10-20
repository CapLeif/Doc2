# Parameter Fields

## Domain and masks

Start from finite field values and valid spatial coordinates. Apply region and dynamical filters independently so a figure title can state exactly what was retained.

| Purpose | Mask |
|---|---|
| Four SD subregions | `combinedArray` in IDs 1–4; see [part.md](part.md) |
| Non-SD/VD regions | IDs 0 and 5 only when the partition data explicitly retains them |
| Large-Ro samples | `abs(Ro) > 1` or `>=1` exactly as specified for that plot |
| Intense negative-eigenvalue samples | `lambda_-/f < -1` |
| Common dual filter | finite base domain AND `abs(Ro)>1` AND `lambda_-/f<-1` |

For a non-SD map, do not use a base mask that has already discarded IDs 0 and 5. For a regional map, combine the chosen region mask with the requested Ro/lambda mask; report the conjunction in the panel label.

## Positive-only quantities

For quantities defined as magnitudes or selected from a positive tail, use a sequential color map and a zero lower bound. In particular, normalized strain magnitude is \(\sigma/|f|\ge0\); M2 frontal-area threshold maps display the nonnegative high-M2 tail. Never use a diverging color scale centered at zero for these positive-only fields.

Set and record the upper color limit per field and output set. Prefer one shared upper limit across panels intended for comparison. When no scientific fixed limit has been established, derive a stated robust limit from finite values (for example, a documented percentile), rather than silently mixing automatic limits. Keep exact physical units in the colorbar label.

Signed variables such as \(\Delta/f\), \(Ro\), and \(\lambda_-/f\) may use a diverging scale centered at zero. State whether normalization uses signed \(f\) or \(|f|\).

## Map layout

- Use x/y in km where the source coordinates are projected distances; preserve equal spatial scaling and the intended domain limits.
- Masked cells remain transparent or visibly neutral; do not paint them as physical zeros. Label the neutral background in the caption if it represents excluded data.
- Use flat cell shading for gridded diagnostics where interpolation would imply unsupported detail. Use a sequential map for positive magnitudes and a diverging map for signed fields.
- Include a colorbar for quantitative maps. Label the field and units; show the exact limits and tick spacing used. Keep a common colorbar scale across directly compared panels.
- If a count is meaningful, include the retained sample count in the title or panel label.

## Source distinctions

`Field_on_1234.m` supplies historical mask and field-specific examples, including a four-region SD mask and `abs(Ro)`/lambda filters. Its old symmetric ranges for strain and M2 are not carried over to positive-only magnitude views. `Filaments/plot_day40_density_ro_lambda_filters.m` supplies an example of finite masks, multi-panel filter comparison, neutral masked background, equal map aspect, and high-resolution white-background export.

