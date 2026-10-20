---
name: mitgcm-llc4320-processing
description: Process compatible local MITgcm LLC4320 surface fields, Kuroshio frontal lines, regional JPDFs, and eigenmode diagnostics using the packaged project algorithms.
---

Use this Skill for the MITgcm LLC4320 data workflows documented in this package.

1. Read [the workflow reference](references/mitgcm-llc4320-processing-flow.md) and inspect the chosen script before running it.
2. Identify the data route first: fixed 960×931×40 Kuroshio binaries, generic 2-D big-endian binaries, or the West Atlantic NetCDF layout. Do not mix their readers or assume shapes/levels.
3. Confirm byte order, Fortran reshape order, timestamp, variable dimensions, face staggering, grid metrics, wet masks, and rotation fields against the actual files.
4. For a new grid or snapshot, validate a small run against a trusted reference before expanding to a batch unless the user explicitly asks to run all.
5. Preserve the implemented formula, mask, plotting range, normalization, and selected region unless the user requests a change. Record real counts and parameters from the current run.
6. Use the included `scripts/` source and install `scripts/requirements.txt`. Model binaries and NetCDF files must be supplied locally; they are not part of this Skill or repository.

