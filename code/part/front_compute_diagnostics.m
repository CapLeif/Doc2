function diagnostics = front_compute_diagnostics(data, faMask)
%FRONT_COMPUTE_DIAGNOSTICS Rotate the velocity gradient into frontal coordinates.

required = {'pupx','pupy','pvpx','pvpy','pbhpx','pbhpy','lonir','latir'};
for k = 1:numel(required)
    if ~isfield(data, required{k})
        error('Front:MissingVariable', 'Required field missing: %s', required{k});
    end
end

ux = data.pupx;
uy = data.pupy;
vx = data.pvpx;
vy = data.pvpy;
gradB = hypot(data.pbhpx, data.pbhpy);
nx = data.pbhpx ./ gradB;
ny = data.pbhpy ./ gradB;
tx = -ny;
ty = nx;

a = nx .* (ux .* nx + uy .* ny) + ny .* (vx .* nx + vy .* ny);
b = nx .* (ux .* tx + uy .* ty) + ny .* (vx .* tx + vy .* ty);
c = tx .* (ux .* nx + uy .* ny) + ty .* (vx .* nx + vy .* ny);
d = tx .* (ux .* tx + uy .* ty) + ty .* (vx .* tx + vy .* ty);

sigma = hypot(a - d, b + c);
valid = faMask & isfinite(a) & isfinite(b) & isfinite(c) & isfinite(d) & sigma > eps;
fields = {a, b, c, d};
normalized = cell(size(fields));
for k = 1:numel(fields)
    normalized{k} = fields{k} ./ sigma;
    normalized{k}(~valid) = NaN;
end
frontogenesis = -a ./ sigma;
H = -4 .* a .* d ./ sigma.^2;
C = 4 .* b .* c ./ sigma.^2;
R = H + C;
frontogenesis(~valid) = NaN;
H(~valid) = NaN;
C(~valid) = NaN;
R(~valid) = NaN;

traceError = (a + d) - (ux + vy);
detError = (a .* d - b .* c) - (ux .* vy - uy .* vx);
strainOriginal = hypot(ux - vy, uy + vx);
strainError = sigma - strainOriginal;
vorticityError = (c - b) - (vx - uy);

diagnostics = struct();
diagnostics.a = a;
diagnostics.b = b;
diagnostics.c = c;
diagnostics.d = d;
diagnostics.a_over_sigma = normalized{1};
diagnostics.b_over_sigma = normalized{2};
diagnostics.c_over_sigma = normalized{3};
diagnostics.d_over_sigma = normalized{4};
diagnostics.frontogenesis = frontogenesis;
diagnostics.H = H;
diagnostics.C = C;
diagnostics.R = R;
diagnostics.sigma = sigma;
diagnostics.nx = nx;
diagnostics.ny = ny;
diagnostics.tx = tx;
diagnostics.ty = ty;
diagnostics.valid_mask = valid;
diagnostics.trace_error = traceError;
diagnostics.det_error = detError;
diagnostics.strain_error = strainError;
diagnostics.vorticity_error = vorticityError;
diagnostics.max_trace_error = max(abs(traceError(valid)), [], 'omitnan');
diagnostics.max_det_error = max(abs(detError(valid)), [], 'omitnan');
diagnostics.max_strain_error = max(abs(strainError(valid)), [], 'omitnan');
diagnostics.max_vorticity_error = max(abs(vorticityError(valid)), [], 'omitnan');
end

