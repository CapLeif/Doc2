# Doc2 research workflows

This repository contains portable instructions, reference notes, and source code for two workflows:

1. the approved Kuroshio Extension Day 40 FL01 combined-figure style and its FL02–FL50 batch application;
2. the MITgcm LLC4320 processing scripts used in this local project.

## Contents

- `docs/day40-fl01-figure-spec.md` — figure layout, visual encodings, processing parameters, and validation facts from the referenced conversation.
- `skills/doc2-day40-fl-figure/` — Codex Skill with its portable reference copy.
- `docs/mitgcm-llc4320-processing-flow.md` — input formats, script boundaries, formulas, and execution sequence.
- `skills/mitgcm-llc4320-processing/` — Codex Skill with its portable reference copy.
- each `skills/.../scripts/` directory — the executable source needed by that Skill, plus its Python requirements.

Raw MITgcm binaries, NetCDF datasets, generated figure collections, and local review documents are intentionally excluded. Supply the data locally and pass its path to the scripts. Do not treat the project-specific defaults as universal LLC4320 conventions.

## Start here

Read the relevant document under `docs/`, then open its matching Skill. On a new computer, copy the complete Skill folder into `%USERPROFILE%\.codex\skills\` (or the configured Codex skills directory). Each folder contains its reference document and required Python source; MITgcm input data must still be supplied separately.

The Day 40 batch that followed the approved FL01 used 10% top-M2 selection, 10 km minimum FA major axis, a strict `>20 km` length gate, at most 50 plotted FLs, 0.5° JPDF bins, a shared 10 km arrow half-length scale, and selected IDs FL02–FL50. The reviewed FL01 was retained. These values describe that run, not every future dataset.

