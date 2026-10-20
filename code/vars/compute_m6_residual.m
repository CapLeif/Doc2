function result = compute_m6_residual(matFile, day)
%COMPUTE_M6_RESIDUAL Compute (sigma^2 - Delta^2 - zeta^2)/f^2 in filaments.

if nargin < 2
    day = 40;
end

if exist(matFile, 'file') ~= 2
    error('M6:MissingInput', 'Input MAT file not found: %s', matFile);
end

requiredVars = { ...
    'Strain_f', 'Diver_f', 'Ro', 'eig2_T_f', 'combinedArray', ...
    'pupx', 'pupy', 'pvpx', 'pvpy', 'f', 'lonir', 'latir'};

S = load(matFile, requiredVars{:});
for k = 1:numel(requiredVars)
    if ~isfield(S, requiredVars{k})
        error('M6:MissingVariable', 'Required variable missing: %s', requiredVars{k});
    end
end

zeta_tensor_f = (S.pvpx - S.pupy) ./ S.f;
theta_all = compute_theta_from_zeta_sigma(zeta_tensor_f, S.Strain_f);
residual_f2 = S.Strain_f.^2 - S.Diver_f.^2 - zeta_tensor_f.^2;
residual_ro_f2 = S.Strain_f.^2 - S.Diver_f.^2 - S.Ro.^2;
detLambda_f2 = (S.pupx .* S.pvpy - S.pvpx .* S.pupy) ./ (S.f.^2);
identity_error = residual_f2 + 4 .* detLambda_f2;

region_mask = ismember(S.combinedArray, 1:4);
finite_main = isfinite(S.Strain_f) & isfinite(S.Diver_f) & ...
    isfinite(S.Ro) & isfinite(S.eig2_T_f) & isfinite(zeta_tensor_f);
filament_mask = region_mask & abs(S.Ro) > 1 & S.eig2_T_f < -1 & finite_main;

det_finite_mask = filament_mask & isfinite(detLambda_f2) & isfinite(identity_error);
if nnz(det_finite_mask) == 0
    identity_max_abs = NaN;
else
    identity_max_abs = max(abs(identity_error(det_finite_mask)));
end

denom = S.Strain_f.^2 + S.Diver_f.^2 + zeta_tensor_f.^2;
relative_abs_residual = abs(residual_f2) ./ denom;
relative_abs_residual(~isfinite(relative_abs_residual) | denom <= eps) = NaN;

summary_table = build_summary_table(residual_f2, detLambda_f2, ...
    relative_abs_residual, identity_error, filament_mask, S.combinedArray);

result = struct();
result.day = day;
result.mat_file = matFile;
result.residual_f2 = residual_f2;
result.residual_ro_f2 = residual_ro_f2;
result.detLambda_f2 = detLambda_f2;
result.identity_error = identity_error;
result.relative_abs_residual = relative_abs_residual;
result.filament_mask = filament_mask;
result.region_mask = region_mask;
result.combinedArray = S.combinedArray;
result.lonir = S.lonir;
result.latir = S.latir;
result.Strain_f = S.Strain_f;
result.Diver_f = S.Diver_f;
result.Ro = S.Ro;
result.zeta_tensor_f = zeta_tensor_f;
result.eig2_T_f = S.eig2_T_f;
result.theta_all = theta_all;
result.summary_table = summary_table;
result.n_filament = nnz(filament_mask);
result.identity_max_abs = identity_max_abs;
end

function T = build_summary_table(residual_f2, detLambda_f2, relative_abs_residual, ...
    identity_error, filament_mask, combinedArray)

group_names = [{'filaments_all'}, arrayfun(@(r) sprintf('region_%d', r), 1:4, 'UniformOutput', false)];
N = zeros(numel(group_names), 1);
mean_residual_f2 = nan(numel(group_names), 1);
median_residual_f2 = nan(numel(group_names), 1);
iqr_residual_f2 = nan(numel(group_names), 1);
rms_residual_f2 = nan(numel(group_names), 1);
p05_residual_f2 = nan(numel(group_names), 1);
p25_residual_f2 = nan(numel(group_names), 1);
p75_residual_f2 = nan(numel(group_names), 1);
p95_residual_f2 = nan(numel(group_names), 1);
mean_abs_residual_f2 = nan(numel(group_names), 1);
median_abs_residual_f2 = nan(numel(group_names), 1);
frac_abs_le_0p1 = nan(numel(group_names), 1);
frac_abs_le_0p25 = nan(numel(group_names), 1);
frac_abs_le_0p5 = nan(numel(group_names), 1);
frac_rel_le_0p05 = nan(numel(group_names), 1);
frac_rel_le_0p10 = nan(numel(group_names), 1);
mean_detLambda_f2 = nan(numel(group_names), 1);
median_detLambda_f2 = nan(numel(group_names), 1);
max_abs_identity_error = nan(numel(group_names), 1);

for i = 1:numel(group_names)
    if i == 1
        mask = filament_mask;
    else
        region_id = i - 1;
        mask = filament_mask & combinedArray == region_id;
    end

    residual_values = residual_f2(mask);
    residual_values = residual_values(isfinite(residual_values));
    abs_residual_values = abs(residual_values);
    rel_values = relative_abs_residual(mask);
    rel_values = rel_values(isfinite(rel_values));
    det_values = detLambda_f2(mask);
    det_values = det_values(isfinite(det_values));
    identity_values = identity_error(mask);
    identity_values = identity_values(isfinite(identity_values));

    N(i) = numel(residual_values);
    if N(i) == 0
        continue;
    end

    pct = local_percentile(residual_values, [5, 25, 50, 75, 95]);
    mean_residual_f2(i) = mean(residual_values);
    median_residual_f2(i) = pct(3);
    iqr_residual_f2(i) = pct(4) - pct(2);
    rms_residual_f2(i) = sqrt(mean(residual_values.^2));
    p05_residual_f2(i) = pct(1);
    p25_residual_f2(i) = pct(2);
    p75_residual_f2(i) = pct(4);
    p95_residual_f2(i) = pct(5);
    mean_abs_residual_f2(i) = mean(abs_residual_values);
    median_abs_residual_f2(i) = median(abs_residual_values);
    frac_abs_le_0p1(i) = mean(abs_residual_values <= 0.1);
    frac_abs_le_0p25(i) = mean(abs_residual_values <= 0.25);
    frac_abs_le_0p5(i) = mean(abs_residual_values <= 0.5);
    frac_rel_le_0p05(i) = mean(rel_values <= 0.05);
    frac_rel_le_0p10(i) = mean(rel_values <= 0.10);

    if ~isempty(det_values)
        mean_detLambda_f2(i) = mean(det_values);
        median_detLambda_f2(i) = median(det_values);
    end
    if ~isempty(identity_values)
        max_abs_identity_error(i) = max(abs(identity_values));
    end
end

T = table(group_names(:), N, mean_residual_f2, median_residual_f2, ...
    iqr_residual_f2, rms_residual_f2, p05_residual_f2, p25_residual_f2, ...
    p75_residual_f2, p95_residual_f2, mean_abs_residual_f2, ...
    median_abs_residual_f2, frac_abs_le_0p1, frac_abs_le_0p25, ...
    frac_abs_le_0p5, frac_rel_le_0p05, frac_rel_le_0p10, ...
    mean_detLambda_f2, median_detLambda_f2, max_abs_identity_error, ...
    'VariableNames', {'group', 'N', 'mean_residual_f2', ...
    'median_residual_f2', 'iqr_residual_f2', 'rms_residual_f2', ...
    'p05_residual_f2', 'p25_residual_f2', 'p75_residual_f2', ...
    'p95_residual_f2', 'mean_abs_residual_f2', ...
    'median_abs_residual_f2', 'frac_abs_le_0p1', 'frac_abs_le_0p25', ...
    'frac_abs_le_0p5', 'frac_rel_le_0p05', 'frac_rel_le_0p10', ...
    'mean_detLambda_f2', 'median_detLambda_f2', 'max_abs_identity_error'});
end

function q = local_percentile(x, p)
x = sort(x(:));
n = numel(x);
q = nan(size(p));
if n == 0
    return;
end
for k = 1:numel(p)
    pos = 1 + (n - 1) * p(k) / 100;
    lo = floor(pos);
    hi = ceil(pos);
    if lo == hi
        q(k) = x(lo);
    else
        q(k) = x(lo) + (pos - lo) * (x(hi) - x(lo));
    end
end
end

