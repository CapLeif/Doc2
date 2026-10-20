# Variables

Only variable-specific definitions are listed here. Shared units and notation are in [README.md](README.md).

## Velocity gradient and strain

For horizontal velocity \(\mathbf{u}=(u,v)\),

\[
\Lambda=\nabla_h\mathbf{u}=\begin{bmatrix}u_x&u_y\\v_x&v_y\end{bmatrix},\qquad
\Lambda^\mathsf{T}=\begin{bmatrix}u_x&v_x\\u_y&v_y\end{bmatrix},\qquad
S=\tfrac12(\Lambda+\Lambda^\mathsf{T}).
\]

The source fields `pupx`, `pupy`, `pvpx`, and `pvpy` represent \(u_x,u_y,v_x,v_y\), respectively. `pbhpx` and `pbhpy` represent \(b_x,b_y\).

| Symbol | Strict definition | Units / notes |
|---|---|---|
| \(\Delta\), `Diver_f` before normalization | \(u_x+v_y=\operatorname{tr}(\Lambda)\) | s\(^{-1}\); divergence |
| \(\zeta\) | \(v_x-u_y\) | s\(^{-1}\); relative vertical vorticity |
| \(\sigma\) | \(\sqrt{(u_x-v_y)^2+(u_y+v_x)^2}\) | s\(^{-1}\); nonnegative strain magnitude |
| \(f\) | Local Coriolis parameter | s\(^{-1}\); may be signed |
| \(Ro\) | \(\zeta/f\) | Dimensionless and signed. Reuse the saved project field for masks so its grid and validity mask stay aligned with the other inputs. |
| \(\Omega\) | \(\sigma^2-\zeta^2\) | s\(^{-2}\); \(\Omega\ge0\) gives real eigenvalues, \(\Omega<0\) gives a complex-conjugate pair |
| \(\lambda_\pm\) | \(\tfrac12(\Delta\pm\sqrt{\Omega})\), when \(\Omega\ge0\) | s\(^{-1}\); algebraically smaller/larger eigenvalues |
| \(\lambda_- / f\) | \(\tfrac12(\Delta/f-\sqrt{\Omega}/|f|)\) on the real branch | Dimensionless and signed; in the lambda-JPDF mapping use the explicit branch coordinates in `jpdf.md` |

`Strain_f` in legacy input files is a stored normalized field and can be sign-sensitive if normalized by signed \(f\). For a nonnegative normalized magnitude use \(\sigma/|f|\), not the ambiguous legacy name.

## Orientation variables

Let \(V_-\), \(V_+\) be the real eigenvectors of \(\Lambda^\mathsf{T}\) ordered by ascending eigenvalue, and \(V_{s-}\), \(V_{s+}\) be eigenvectors of \(S\) selected by negative/positive signed eigenvalue.

For unoriented axes, define \(\phi(a,b)=\cos^{-1}(|\hat a\cdot\hat b|)\), with normalized vectors and range \([0^\circ,90^\circ]\).

| Variable | Definition | Validity |
|---|---|---|
| \(\theta\) | Directed angle between the real eigendirections of \(\Lambda^\mathsf{T}\), equivalently \(\cos^{-1}(-\zeta/\sigma)\) | Directed angle in \([0^\circ,180^\circ]\); requires \(\sigma>0\), real eigenstructure, and a clamped cosine within numerical tolerance |
| \(\alpha\) | \(\phi(V_-,V_{s-})\) | Requires a real, distinct eigenpair of \(\Lambda^\mathsf{T}\) and one negative/one positive eigenvalue of \(S\) |
| \(\gamma\) | \(\phi(V_{s-},\nabla_h b)\) | Requires nonzero finite buoyancy gradient and valid strain compression axis |
| \(\epsilon\) | \(\phi(V_-,\nabla_h b)\) | Requires nonzero finite buoyancy gradient and valid negative mode of \(\Lambda^\mathsf{T}\) |

The geometry check is \(\min(\theta,180^\circ-\theta)+2\alpha=90^\circ\), up to numerical error. Where the paired historical notation calls the same angle \(\beta\), \(\alpha=\beta\); this project keeps only \(\alpha\) as the variable name.

Signed variants `gamma_signed` and `epsilon_signed`, where used, retain signed acute-axis angles; their sign is determined from the local anticyclonic/cyclonic sense and local \(f\). Do not substitute them for unsigned \(\gamma\) or \(\epsilon\).

## Buoyancy-gradient forcing

With \(b_x=\partial_x b\), \(b_y=\partial_y b\),

\[
F_{\mathrm{hor},b}=-u_{k,j}b_kb_j,
\qquad T_b=\frac{F_{\mathrm{hor},b}}{|\nabla_h b|^2}.
\]

The implementation expands this as

\[
F_{\mathrm{hor},b}=-\left(u_xb_x^2+(u_y+v_x)b_xb_y+v_yb_y^2\right),
\qquad T_b=F_{\mathrm{hor},b}/(b_x^2+b_y^2).
\]

Require finite gradients and \(|\nabla_hb|^2>0\); `Tb/f` is a signed normalized rate.

## M2

`M2` is the supplied project field used for frontal-area thresholding and as the background field in FL figures. The extraction sources do not define its upstream physical formula, so this convention does not invent one. Record upstream units and derivation with the data product that creates `M2`. In the extraction workflow, the accepted high-M2 tail is treated as nonnegative.

## Strain determinant identity

For \(\det(\Lambda)=u_xv_y-u_yv_x\),

\[
\sigma^2-\Delta^2-\zeta^2=-4\det(\Lambda),
\qquad
\frac{\zeta^2+\Delta^2}{\sigma^2}-\frac{4\det(\Lambda)}{\sigma^2}-1=0.
\]

Use the same velocity-gradient tensor for all terms in a numerical check.

