# VGL Partitions

The six partition labels are defined by the project as follows. V, G, and L mean **Vortical**, **Genesis**, and **Lysis**. The numeric ID-to-label mapping below is authoritative for this project; check that a saved partition field uses the same encoding before applying it.

## Coordinates and branch

Use \(x=\Delta/f\) and \(\Omega_f=(\sigma/f)^2-(\zeta/f)^2\) for the algebraic tests. On the real branch \(\Omega_f\ge0\), \(q=\sqrt{\Omega_f}\). The lambda-JPDF display uses a signed vertical coordinate; its exact branch convention is in [jpdf.md](jpdf.md). The physical inequalities below are based on \(\Delta\), \(\sqrt{\Omega}\), and whether strain or vorticity dominates.

## Letter meanings and order

- **V** marks a vorticity-dominated regime.
- For real eigenvalues, **G** denotes Genesis and corresponds to the negative eigenvalue; **L** denotes Lysis and corresponds to the positive eigenvalue.
- When a real-eigenvalue region has one negative and one positive eigenvalue, the first letter is the process associated with the eigenvalue of larger absolute magnitude. Thus **LG** means the positive/Lysis eigenvalue is stronger, while **GL** means the negative/Genesis eigenvalue is stronger.
- In the two vorticity-dominated regions, V comes first. The second letter distinguishes the divergent/Lysis branch (\(\Delta>0\), L) from the convergent/Genesis branch (\(\Delta<0\), G).
- Equal-magnitude eigenvalues and equality surfaces are boundaries; do not assign them by an arbitrary tie-break.

## Regime table

| ID | Label | Meaning | Flow regime and criteria |
|---:|:---:|---|---|
| 0 | **VL** | Vortical + Lysis | Divergent vorticity-dominated complex eigenpair; \(\Delta>0\), \(\zeta^2>\sigma^2\) |
| 1 | **LL** | Lysis + Lysis | Divergence-dominated strained flow; \(\Delta>0\), \(\sigma^2>\zeta^2\), \(|\Delta|>\sqrt{\Omega}\) |
| 2 | **LG** | Lysis stronger than Genesis | Strain-dominated divergent flow; \(\Delta>0\), \(\sigma^2>\zeta^2\), \(\sqrt{\Omega}>|\Delta|\) |
| 3 | **GL** | Genesis stronger than Lysis | Strain-dominated convergent flow; \(\Delta<0\), \(\sigma^2>\zeta^2\), \(\sqrt{\Omega}>|\Delta|\) |
| 4 | **GG** | Genesis + Genesis | Convergence-dominated strained flow; \(\Delta<0\), \(\sigma^2>\zeta^2\), \(|\Delta|>\sqrt{\Omega}\) |
| 5 | **VG** | Vortical + Genesis | Convergent vorticity-dominated complex eigenpair; \(\Delta<0\), \(\zeta^2>\sigma^2\) |

Regions 1–4 are strain dominated (SD); regions 0 and 5 are vorticity dominated (VD). Equality boundaries \(\Delta=0\), \(\sigma^2=\zeta^2\), and \(|\Delta|=\sqrt{\Omega}\) are boundaries, not extra classes.

Some saved combinedArray fields encode only the four real-eigenvalue SD subregions (IDs 1–4). Confirm how each saved array was constructed before applying the complete 0–5 mapping; do not infer that absent IDs or other values mean VD.

## Mask semantics

- SD mask: IDs 1–4, or their equivalent real-eigenvalue criteria.
- VD/non-SD mask: IDs 0 and 5 only when the saved partition variable retains the complete six-region encoding.
- A commonly used dual filter for intense frontogenetic samples is finite fields plus \(|Ro|>1\) and \(\lambda_-/f<-1\). These strict inequalities are part of that filter; some field maps use \(|Ro|\ge1\), so retain the threshold specified by the target figure.

## Sources

Project manuscript Table 1 and its accompanying partition discussion; Field_on_1234.m; M6_sigma_delta_zeta_day40/src/front_compute_fl_lambda_jpdf.m; M6_sigma_delta_zeta_day40/src/compute_m6_residual.m.
