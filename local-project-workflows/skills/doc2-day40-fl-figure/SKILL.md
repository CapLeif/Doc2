---
name: doc2-day40-fl-figure
description: Reproduce the approved Kuroshio Extension Day 40 FL combined figure with M2, eigendirection arrows, signed angle/rate diagnostics, and strain-identity panel.
---

Use this Skill when working on the Day 40 FL01–FL50 figure family or a compatible Kuroshio LLC4320 snapshot.

This package covers the Python Kuroshio LLC4320 figure family. Keep it distinct from the separate MATLAB `εγθ_NewFL` plotting contract in the repository's Doc2 core reference.

1. Read [the figure specification](references/day40-fl01-figure-spec.md) before changing the plot or running a batch.
2. Use the included repository source `scripts/kuroshio_fl_theta_gamma_epsilon.py` as the executable calculation and drawing implementation. Its `plot_preview`/`draw_mode` functions define the figure; its `run_batch` path defines the inventory and export behavior.
3. Confirm the input snapshot, expected binary layout, and output directory. Do not infer matching FL IDs from an earlier timestamp.
4. Preserve the approved three-panel layout, colors, 10 km shared arrow scale, direct square-root mapping, display-only centerline smoothing, mapped λ/rate bands, and identity residual gate unless the user requests a change.
5. For a reviewed-preview workflow, render the requested preview first and wait for the user's review before the remaining batch. If the user explicitly asks to render all directly, proceed.
6. Keep raw model binaries out of the repository. Report the actual IDs, number of files written, and any skipped/non-drawable lines.

