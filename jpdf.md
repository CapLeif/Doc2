# Lambda-JPDF

This file gives the full-field lambda-JPDF contract. Plot-specific ranges can differ for conditional or cropped products and must be stated with that output.

## Coordinates

Let

\[
x=\Delta/f,
\qquad \Omega_f=(\sigma/|f|)^2-(\zeta/f)^2.
\]

Display the vertical coordinate as

\[
y=\begin{cases}
-\sqrt{\Omega_f}, & \Omega_f\ge0 \quad(\text{real eigenvalues}),\\
+\sqrt{-\Omega_f}, & \Omega_f<0 \quad(\text{complex pair}).
\end{cases}
\]

Thus the real branch lies below zero and the complex branch above zero. The plotted eigenvalue contours use the real branch only:

\[
\lambda_-/f=\tfrac12(x+y),\qquad
\lambda_+/f=\tfrac12(x-y).
\]

## Histogram and display

- Build uniform two-dimensional bins and normalize counts as \(p_{ij}=n_{ij}/(N_{\mathrm{in}}\,\Delta x\Delta y)\), where \(N_{\mathrm{in}}\) is the number of valid samples inside the plotted range. This makes the displayed in-range density integrate to one.
- The full-field reference script uses \([-10,10]\) on each axis with bin width 0.05. Treat these as that figure's configuration, not a universal requirement.
- Map positive density logarithmically with \(\log_{10}p\); zero-count bins are white. The reference uses a low reference level of 0.001 and labels 0.001, 0.01, 0.1, and 1; the upper color limit is the observed maximum.
- Preserve equal data scaling (`axis equal`), white background, Arial axes, outward ticks, visible zero axes, and 600 dpi export for this reference style.
- Draw dashed \(\lambda_-/f\) contours at integer levels -10 to 0 in the real branch. Draw \(\lambda_+/f\) contours at the same levels only in the specified third-quadrant subregion. Keep contour labels legible and distinguish the two families by color.
- Label axes with the exact signed normalizations: \(\Delta/f\) and the branch-wise \(-\sqrt{\Omega_f}\) / \(+\sqrt{-\Omega_f}\). Do not shorten the latter to an ambiguous `sqrtOmega` label.

## Part labels

Use the names and boundaries in [part.md](part.md): iDv+, Ds+, Sd−, Sc−, Cs−, and iCv−. Boundary points should not receive a class unless the selected code explicitly defines a tie rule. The strict interior labels must remain consistent with the signs of \(\Delta\) and the real/complex branch.

## Conditional products

For conditional statistics, store and report counts per bin, bin width, minimum count, sample mask, and whether the statistic is a mean or median. Bins below the declared minimum count are missing, not zero. Existing sigma-zeta conditional scripts use 0.05 bins and minimum count 5; do not silently transfer those settings to a different lambda-JPDF output.

## Source

`M6_sigma_delta_zeta_day40/src/compute_frontogenesis_diagnostics.m`, `M6_sigma_delta_zeta_day40/scripts/run_fullfield_lambda_jpdf_day40.m`, and the project manuscript's lambda-JPDF figure and partition discussion.

