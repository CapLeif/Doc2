# Front Extraction and Arrow Figures

## Candidate frontal areas

1. Restrict to the selected domain and finite M2 values.
2. Select the requested upper M2 percentile and form a threshold mask.
3. Find 8-connected components. For a region workflow retain the selected component(s); for full-field discovery apply the declared minimum major-axis length and preserve an audit of rejected candidates.
4. State percentile, spatial resolution, connectivity, size gates, and whether the operation is ROI or full-field. These settings vary by workflow and are not universal constants.

## Constrained ridge extraction

Use the M2 field with buoyancy-gradient components \(b_x,b_y\). Smooth M2, compute its gradient and Hessian, and define unit buoyancy-gradient normal \(n=\nabla_hb/|\nabla_hb|\). Candidate ridge points satisfy negative normal curvature and a small normal-offset correction; keep candidates inside the accepted frontal-area mask. Thin the support, split branch points, connect neighboring points only when the path aligns with the local tangent, and retain continuous paths above the declared minimum length. Resample along arc length and recheck mask membership and negative normal curvature.

The FL01/NewFL generating workflow uses top-15% M2 connected regions, 10 km minimum component major axis, 10 km minimum final FL length, 45° maximum tangent mismatch, and Gaussian smoothing with grid sigma 1 for ridge extraction. Its FL01 contract must not be confused with the separate full-field top-10% workflow.

## FL01 arrow-panel contract

- Background: M2 on the native x/y grid, `summer` sequential palette, labeled colorbar, equal spatial aspect, white figure background.
- Overlay: extracted FL as a black line; show each eigen-sample location as a small white-filled black marker.
- Matrix roles: `Vs−`, `Vs+` are signed compression/extension eigenvectors of \(S\); `V−`, `V+` are negative/positive real-eigenvalue eigenvectors of \(\Lambda^\mathsf{T}\). Select by signed eigenvalue, not magnitude. Keep eigenvector signs continuous between neighboring samples when rendering a sequence.
- Legend colors in the FL01 generating script: `Vs−` gray, `Vs+` dark gray, `V−` blue, `V+` red. Negative modes are drawn inward toward the sample point as opposing arrows; positive modes point outward from it. Disable automatic quiver scaling.
- Scale arrow half-length from eigenvalue magnitude relative to one global maximum across the compared samples: \(L=L_{max}(|\lambda|/\lambda_{max})^{1/p}\). The FL01 source uses a global scale, a 10 km nominal maximum half-length, and square-root scaling; do not replace this with per-front scaling when comparing fronts.
- Resample arrows at roughly 3 km spacing while retaining both endpoints. The NewFL combined export uses a three-row/two-column landscape layout at 600 dpi: map at left, along-FL angle/diagnostic panel at upper right, and strain-identity panel at lower right.
- For display, smooth the drawn FL at 0.25 km intervals with a 2 km Gaussian scale while fixing both endpoints; retain the unsmoothed extracted ridge for calculations.
- Use one global eigenvalue scale across all plotted fronts. The reference uses a 10 km maximum half-length and square-root scaling (power exponent 2), with quiver autoscaling disabled. Negative modes point inward; positive modes point outward.
- The upper-right panel uses distance from the upper FL endpoint and plots θ with signed γ and signed ε. Its angle axis spans −90° to 180°; the normalized eigenvalue and frontogenetic-rate bands map to 60°–120° on that panel.
- The lower-right panel shows the strain identity as stacked contributions: -4det(Λ)/σ², then ζ²/σ², then δ²/σ². The sum is 1 within numerical tolerance; expand the vertical limits to show negative determinant compensation when present.

## Sampling and display safeguards

The referenced FL01 script samples its own input inventory; other front workflows use different sample fractions. Record the fractions, any relocation radius, endpoint handling, scale exponent, and max arrow length with the figure. Skip invalid/complex eigenpairs instead of drawing a fabricated direction. Include the legend mapping and ensure the crop leaves room around the longest arrows.

## Sources

The reference figure is `FL_lambdaJPDF/εγθ_NewFL/combined_figures/FL01_arrows_theta_gamma_epsilon.png`. Its generating code is `M6_sigma_delta_zeta_day40/scripts/run_day40_newfl_theta_gamma_epsilon.m`; shared extraction logic is in `M6_sigma_delta_zeta_day40/src/front_extract_fullfield_fas.m` and `front_extract_ridges.m`.

