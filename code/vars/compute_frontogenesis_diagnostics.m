function diagnostics = compute_frontogenesis_diagnostics(data)
%COMPUTE_FRONTOGENESIS_DIAGNOSTICS Compute horizontal buoyancy frontogenesis.
%
% F_hor,b = -u_{k,j} b_{,k} b_{,j} and
% T_b = F_hor,b / |grad_h b|^2 follow Barkan et al. (2019), Eqs. (1),
% (13a), and (15c), and Yu et al. (2024), Methods Eq. (1).

required = {'pupx','pupy','pvpx','pvpy','pbhpx','pbhpy','f'};
for k = 1:numel(required)
    if ~isfield(data, required{k})
        error('Frontogenesis:MissingVariable', ...
            'Required field missing: %s', required{k});
    end
end

ux = data.pupx;
uy = data.pupy;
vx = data.pvpx;
vy = data.pvpy;
bx = data.pbhpx;
by = data.pbhpy;
f = data.f;

referenceSize = size(ux);
fields = {uy, vx, vy, bx, by, f};
for k = 1:numel(fields)
    if ~isequal(size(fields{k}), referenceSize)
        error('Frontogenesis:SizeMismatch', ...
            'All velocity-gradient, buoyancy-gradient, and f fields must have the same size.');
    end
end

grad_b_squared = bx.^2 + by.^2;
q_x = -(ux .* bx + vx .* by);
q_y = -(uy .* bx + vy .* by);
frontogenesis_function = q_x .* bx + q_y .* by;
frontogenetic_rate = frontogenesis_function ./ grad_b_squared;
frontogenetic_rate_f = frontogenetic_rate ./ f;

delta_f = (ux + vy) ./ f;
sigma_f = hypot(ux - vy, uy + vx) ./ abs(f);
zeta_f = (vx - uy) ./ f;
omega_f2 = sigma_f.^2 - zeta_f.^2;

lambda_minus_f = nan(referenceSize);
realEigenMask = omega_f2 >= 0 & isfinite(omega_f2) & isfinite(delta_f);
lambda_minus_f(realEigenMask) = 0.5 .* (...
    delta_f(realEigenMask) - sqrt(omega_f2(realEigenMask)));

lambda_jpdf_y = nan(referenceSize);
lambda_jpdf_y(realEigenMask) = -sqrt(omega_f2(realEigenMask));
complexEigenMask = omega_f2 < 0 & isfinite(omega_f2) & isfinite(delta_f);
lambda_jpdf_y(complexEigenMask) = sqrt(-omega_f2(complexEigenMask));

valid = isfinite(frontogenesis_function) & isfinite(frontogenetic_rate_f) & ...
    isfinite(delta_f) & isfinite(lambda_jpdf_y) & isfinite(f) & ...
    grad_b_squared > 0 & f ~= 0;

frontogenesis_function(~valid) = NaN;
frontogenetic_rate(~valid) = NaN;
frontogenetic_rate_f(~valid) = NaN;
q_x(~valid) = NaN;
q_y(~valid) = NaN;

diagnostics = struct();
diagnostics.grad_b_squared = grad_b_squared;
diagnostics.q_x = q_x;
diagnostics.q_y = q_y;
diagnostics.frontogenesis_function = frontogenesis_function;
diagnostics.frontogenetic_rate = frontogenetic_rate;
diagnostics.frontogenetic_rate_f = frontogenetic_rate_f;
diagnostics.delta_f = delta_f;
diagnostics.sigma_f = sigma_f;
diagnostics.zeta_f = zeta_f;
diagnostics.omega_f2 = omega_f2;
diagnostics.lambda_minus_f = lambda_minus_f;
diagnostics.lambda_jpdf_y = lambda_jpdf_y;
diagnostics.real_eigen_mask = realEigenMask;
diagnostics.complex_eigen_mask = complexEigenMask;
diagnostics.valid_mask = valid;
end

