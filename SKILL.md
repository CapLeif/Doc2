---
name: doc2-figure-reference
description: Apply Doc2 scientific variables, diagnostics, field masks, lambda-JPDF, and frontal-line extraction/figure conventions when creating or revising related analyses and plots.
---

# Doc2 figure and diagnostics reference

Use this repository as the canonical project reference for Doc2-related scientific diagnostics and figures. Before editing or creating a figure, read the relevant topic documents and inspect the corresponding MATLAB implementation under `code/`.

## Route by task

- Variable definitions and identities: `vars.md`, then `code/vars/`.
- VGL names, criteria, or flow-region filtering: `part.md`, then `code/part/`. Do not infer unconfirmed numeric IDs.
- λ-JPDF: `jpdf.md`, then `code/jpdf/`.
- Parameter maps and filtering: `field.md`, then `code/field/`. A base/domain/non-SD mask must come from the authoritative workflow; `field_filter_mask` applies only finite, Ro, and λ₋/f filters.
- Front extraction or eigenvector arrows: `front.md`, then `code/front/`.

Read `README.md` for shared symbols, units, normalization, sources, applicability, and repository scope. Topic-specific thresholds and styles are not universal unless their documents explicitly say so.

## Use the code responsibly

Treat MATLAB files as reference implementations with caller-supplied arrays and project data structures. Read each function's signature, required fields, and validity handling before reuse. Keep machine paths and dataset-specific orchestration in the consuming project; do not add data or generated figure batches here.

When a requested plot is not covered by these references, identify the nearest documented convention and state any new choice clearly instead of presenting it as an established Doc2 rule. Update the relevant document and code when the project adopts a new shared convention.



## MITgcm LLC4320 companion workflows

For the Python Kuroshio FL01 figure and its batch convention, read `local-project-workflows/docs/day40-fl01-figure-spec.md` and use `local-project-workflows/skills/doc2-day40-fl-figure/`.

For the broader MITgcm LLC4320 processing flow, read `local-project-workflows/docs/mitgcm-llc4320-processing-flow.md` and use `local-project-workflows/skills/mitgcm-llc4320-processing/`.

These Skills are separate from the MATLAB Doc2 workflows above. Each includes its supporting reference and Python source; model datasets remain external.