function theta = compute_theta_from_zeta_sigma(zeta, sigma)
%COMPUTE_THETA_FROM_ZETA_SIGMA Compute theta from cos(theta) = -zeta/sigma.

if ~isequal(size(zeta), size(sigma))
    error('M6:ThetaSizeMismatch', 'zeta and sigma must have identical sizes.');
end

theta = nan(size(zeta));
valid = isfinite(zeta) & isfinite(sigma) & abs(sigma) > eps;
cos_theta = -zeta(valid) ./ sigma(valid);
cos_theta = max(-1, min(1, cos_theta));
theta(valid) = acosd(cos_theta);
end

