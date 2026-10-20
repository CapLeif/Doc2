function angles = front_compute_fl_angles(ridge, data)
%FRONT_COMPUTE_FL_ANGLES Compute four orientation angles at every FL point.
%
% theta: directed A' eigenvector angle from cos(theta) = -zeta/sigma.
% alpha: principal angle between blue (lambda_minus of A') and S compression.
% beta:  principal angle between red (lambda_plus of A') and S extension.
% gamma: principal angle between Vs- (S compression) and grad_h(b).
% epsilon: principal angle between V- (negative A' mode) and grad_h(b).
% gamma_signed/epsilon_signed: the same acute axis angles with positive
% values in the local anticyclonic sense and negative values in the local
% cyclonic sense. The sign uses the local Coriolis parameter.
%
% The S compression/extension directions are selected by signed eigenvalue,
% never by eigenvalue magnitude. Principal angles use abs(dot), because an
% eigenvector is an unoriented axis: v and -v represent the same direction.

required = {'lonir','latir','pupx','pupy','pvpx','pvpy','pbhpx','pbhpy'};
for k = 1:numel(required)
    if ~isfield(data, required{k})
        error('Front:MissingVariable', 'Required field missing: %s', required{k});
    end
end

% Accept either stored orientation, but always diagnose from the upper end
% toward the lower end. New curated inventories are already top-to-bottom.
x = ridge.x(:);
y = ridge.y(:);
if y(1)<y(end) || (y(1)==y(end) && x(1)>x(end))
    x=flip(x); y=flip(y);
end
arcFromTop=[0;cumsum(hypot(diff(x),diff(y)))];
n = numel(x);

ux = interp2(data.lonir, data.latir, data.pupx, x, y, 'linear');
uy = interp2(data.lonir, data.latir, data.pupy, x, y, 'linear');
vx = interp2(data.lonir, data.latir, data.pvpx, x, y, 'linear');
vy = interp2(data.lonir, data.latir, data.pvpy, x, y, 'linear');
bx = interp2(data.lonir, data.latir, data.pbhpx, x, y, 'linear');
by = interp2(data.lonir, data.latir, data.pbhpy, x, y, 'linear');
if isfield(data,'f')
    coriolisF = interp2(data.lonir,data.latir,data.f,x,y,'linear');
else
    coriolisF = nan(n,1);
end

theta = nan(n, 1);
alpha = nan(n, 1);
beta = nan(n, 1);
gamma = nan(n, 1);
epsilon = nan(n, 1);
gammaSigned = nan(n, 1);
epsilonSigned = nan(n, 1);
sigma = hypot(ux - vy, uy + vx);
zeta = vx - uy;
divergence = ux + vy;
detLambda = ux .* vy - uy .* vx;
sigmaSquared = sigma .^ 2;
zeta2OverSigma2 = zeta .^ 2 ./ sigmaSquared;
divergence2OverSigma2 = divergence .^ 2 ./ sigmaSquared;
fourDetLambdaOverSigma2 = 4 .* detLambda ./ sigmaSquared;
strainIdentityResidual = zeta2OverSigma2 + divergence2OverSigma2 - ...
    fourDetLambdaOverSigma2 - 1;
gradBSquared = bx.^2 + by.^2;
frontogeneticRate = nan(n,1);
validRate = isfinite(ux) & isfinite(uy) & isfinite(vx) & isfinite(vy) & ...
    isfinite(bx) & isfinite(by) & gradBSquared > eps;
frontogeneticRate(validRate) = -( ...
    ux(validRate).*bx(validRate).^2 + ...
    (uy(validRate)+vx(validRate)).*bx(validRate).*by(validRate) + ...
    vy(validRate).*by(validRate).^2) ./ gradBSquared(validRate);
lambdaSCompression = nan(n, 1);
lambdaSExtension = nan(n, 1);
lambdaATMinus = nan(n, 1);
lambdaATPlus = nan(n, 1);
thetaAxisResidual = nan(n, 1);
thetaAlphaIdentityResidual = nan(n, 1);
valid = false(n, 1);

for k = 1:n
    values = [ux(k) uy(k) vx(k) vy(k) bx(k) by(k)];
    if any(~isfinite(values)) || sigma(k) <= eps || hypot(bx(k), by(k)) <= eps
        continue;
    end

    cosTheta = -zeta(k) / sigma(k);
    if abs(cosTheta) > 1 + 1e-10
        % A' has a complex-conjugate eigenpair in this regime.
        continue;
    end
    cosTheta = max(-1, min(1, cosTheta));

    A = [ux(k) uy(k); vx(k) vy(k)];
    [VAT, DAT] = eig(A.');
    evalAT = diag(DAT);
    if max(abs(imag(evalAT))) > 1e-10
        continue;
    end
    [evalAT, orderAT] = sort(real(evalAT), 'ascend');
    VAT = real(VAT(:, orderAT));
    if abs(evalAT(2) - evalAT(1)) <= 100 * eps(max(1, norm(A, 'fro')))
        continue;
    end
    % Theta and the ordered A' eigenvalues remain valid independently of
    % whether S has a signed compression-extension pair.
    theta(k) = acosd(cosTheta);
    lambdaATMinus(k) = evalAT(1);
    lambdaATPlus(k) = evalAT(2);

    S = 0.5 .* (A + A.');
    [VS, DS] = eig(S);
    [evalS, orderS] = sort(real(diag(DS)), 'ascend');
    VS = real(VS(:, orderS));
    if ~(evalS(1) < 0 && evalS(2) > 0)
        % A unique negative compression axis and positive extension axis
        % are required for alpha/beta/gamma.
        continue;
    end

    blue = normalize_vector(VAT(:, 1));
    red = normalize_vector(VAT(:, 2));
    compression = normalize_vector(VS(:, 1));
    extension = normalize_vector(VS(:, 2));
    gradB = normalize_vector([bx(k); by(k)]);

    alpha(k) = principal_axis_angle(blue, compression);
    beta(k) = principal_axis_angle(red, extension);
    gamma(k) = principal_axis_angle(compression, gradB);
    epsilon(k) = principal_axis_angle(blue, gradB);
    if isfinite(coriolisF(k)) && abs(coriolisF(k)) > eps
        gammaSigned(k) = signed_axis_angle_anticyclonic_positive( ...
            gradB, compression, sign(coriolisF(k)));
        epsilonSigned(k) = signed_axis_angle_anticyclonic_positive( ...
            gradB, blue, sign(coriolisF(k)));
    end
    lambdaSCompression(k) = evalS(1);
    lambdaSExtension(k) = evalS(2);

    % Independent geometry check after removing eigenvector sign ambiguity.
    thetaAxis = principal_axis_angle(red, blue);
    thetaAxisResidual(k) = thetaAxis - min(theta(k), 180 - theta(k));
    thetaAlphaIdentityResidual(k) = min(theta(k),180-theta(k)) + 2*alpha(k) - 90;
    valid(k) = true;
end

angles = struct();
angles.grid_index_top_to_bottom = (1:n).';
angles.arc_from_top_km = arcFromTop;
angles.x = x;
angles.y = y;
angles.theta_deg = theta;
angles.alpha_deg = alpha;
angles.beta_deg = beta;
angles.gamma_deg = gamma;
angles.epsilon_deg = epsilon;
angles.gamma_signed_deg = gammaSigned;
angles.epsilon_signed_deg = epsilonSigned;
angles.sigma = sigma;
angles.zeta = zeta;
angles.divergence = divergence;
angles.det_lambda = detLambda;
angles.zeta2_over_sigma2 = zeta2OverSigma2;
angles.divergence2_over_sigma2 = divergence2OverSigma2;
angles.four_det_lambda_over_sigma2 = fourDetLambdaOverSigma2;
angles.strain_identity_residual = strainIdentityResidual;
angles.frontogenetic_rate = frontogeneticRate;
angles.coriolis_f = coriolisF;
angles.lambda_AT_minus_over_f = lambdaATMinus ./ coriolisF;
angles.lambda_AT_plus_over_f = lambdaATPlus ./ coriolisF;
angles.frontogenetic_rate_over_f = frontogeneticRate ./ coriolisF;
angles.lambda_S_compression = lambdaSCompression;
angles.lambda_S_extension = lambdaSExtension;
angles.lambda_AT_minus = lambdaATMinus;
angles.lambda_AT_plus = lambdaATPlus;
angles.theta_axis_residual_deg = thetaAxisResidual;
angles.theta_acute_plus_2alpha_minus_90_deg = thetaAlphaIdentityResidual;
angles.valid = valid;
angles.n_total = n;
angles.n_valid = nnz(valid);
angles.max_abs_theta_axis_residual_deg = ...
    max(abs(thetaAxisResidual(valid)), [], 'omitnan');
angles.max_abs_theta_alpha_identity_residual_deg = ...
    max(abs(thetaAlphaIdentityResidual(valid)), [], 'omitnan');
end

function value = principal_axis_angle(a, b)
value = acosd(max(0, min(1, abs(dot(a, b)))));
end

function value = signed_axis_angle_anticyclonic_positive(reference, axisVector, coriolisSign)
% Fold the directed reference-to-axis rotation into the shortest axis angle.
% This makes the result invariant to the arbitrary eigenvector sign.
ccwAngle = atan2d(reference(1)*axisVector(2)-reference(2)*axisVector(1), ...
    dot(reference,axisVector));
ccwAcute = mod(ccwAngle+90,180)-90;
value = -coriolisSign .* ccwAcute;
end

function v = normalize_vector(v)
v = v ./ norm(v);
end

