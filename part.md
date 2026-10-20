# VGL Partitions

The labels below are the six flow regimes in the project's eigenvalue partition. The manuscript names these regimes and gives their inequalities; `combinedArray` conventionally stores the four real-eigenvalue SD subregions as IDs 1–4. “VGL” is retained only as the requested compact index; no expansion of its letters is inferred.

## Coordinates and branch

Use \(x=\Delta/f\) and \(\Omega_f=(\sigma/f)^2-(\zeta/f)^2\) for the algebraic tests. On the real branch \(\Omega_f\ge0\), \(q=\sqrt{\Omega_f}\). The lambda-JPDF display uses a signed vertical coordinate; its exact branch convention is in [jpdf.md](jpdf.md). The labels' physical inequalities are based on \(\Delta\), \(\sqrt{\Omega}\), and whether strain or vorticity dominates.

## Regime table

| ID | Name | Flow regime | Conditions |
|---:|---|---|---|
| 0 | iDv+ | Divergent vorticity dominated; complex eigenpair | \(\Delta>0\), \(\zeta^2>\sigma^2\) |
| 1 | Ds+ | Divergence-dominated strained flow | \(\Delta>0\), \(\sigma^2>\zeta^2\), \(|\Delta|>\sqrt{\Omega}\) |
| 2 | Sd− | Strain-dominated divergent flow | \(\Delta>0\), \(\sigma^2>\zeta^2\), \(\sqrt{\Omega}>|\Delta|\) |
| 3 | Sc− | Strain-dominated convergent flow | \(\Delta<0\), \(\sigma^2>\zeta^2\), \(\sqrt{\Omega}>|\Delta|\) |
| 4 | Cs− | Convergence-dominated strained flow | \(\Delta<0\), \(\sigma^2>\zeta^2\), \(|\Delta|>\sqrt{\Omega}\) |
| 5 | iCv− | Convergent vorticity dominated; complex eigenpair | \(\Delta<0\), \(\zeta^2>\sigma^2\) |

The manuscript identifies regions 1–4 as SD and regions 0 and 5 as vorticity-dominated (VD). Equality boundaries \(\Delta=0\), \(\sigma^2=\zeta^2\), and \(|\Delta|=\sqrt{\Omega}\) should be treated as boundaries, not assigned by guessed precedence. The four SD IDs above match `combinedArray` use in the existing scripts; verify source mask construction before applying those numeric IDs to a different data product.

## Mask semantics

- SD mask: IDs 1–4, or the equivalent real-eigenvalue criteria with the two dominance subdivisions.
- VD/non-SD mask: IDs 0 and 5 when the saved partition variable retains all six regimes. Do not equate `~ismember(combinedArray,1:4)` with VD unless the array actually encodes IDs 0 and 5 and other values are handled.
- Dual filter often used for intense frontogenetic samples: finite fields plus `abs(Ro)>1` and `lambda_-/f<-1`. The strict inequalities are part of that filter; some older field maps use `abs(Ro)>=1`, so copy the exact requested threshold.

## Sources

Project manuscript Table 1 and its accompanying partition discussion; `Field_on_1234.m`; `M6_sigma_delta_zeta_day40/src/front_compute_fl_lambda_jpdf.m`; `M6_sigma_delta_zeta_day40/src/compute_m6_residual.m`.

