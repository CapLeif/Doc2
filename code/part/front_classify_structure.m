function classification = front_classify_structure(diagnostics, x, y, nBootstrap, blockKm, seed)
%FRONT_CLASSIFY_STRUCTURE Classify hyperbolic versus shear-coupled contributions.

if nargin < 4 || isempty(nBootstrap), nBootstrap = 1000; end
if nargin < 5 || isempty(blockKm), blockKm = 2; end
if nargin < 6 || isempty(seed), seed = 40; end

valid = diagnostics.valid_mask & isfinite(diagnostics.a) & isfinite(diagnostics.b) & ...
    isfinite(diagnostics.c) & isfinite(diagnostics.d);
adAbs = abs(diagnostics.a(valid) .* diagnostics.d(valid));
bcAbs = abs(diagnostics.b(valid) .* diagnostics.c(valid));
denominator = sum(adAbs + bcAbs);
if denominator <= eps || isempty(adAbs)
    error('Front:ClassificationEmpty', 'No finite nonzero FA contributions are available.');
end
FH = sum(adAbs) / denominator;
classicFraction = mean(diagnostics.a(valid) < 0 & diagnostics.d(valid) > 0);
shearFraction = mean(diagnostics.b(valid) .* diagnostics.c(valid) > 0);

xv = x(valid);
yv = y(valid);
blockX = floor((xv - min(xv)) / blockKm);
blockY = floor((yv - min(yv)) / blockKm);
[~, ~, blockIndex] = unique([blockX(:), blockY(:)], 'rows', 'stable');
nBlocks = max(blockIndex);
blockAD = accumarray(blockIndex, adAbs, [nBlocks 1], @sum, 0);
blockBC = accumarray(blockIndex, bcAbs, [nBlocks 1], @sum, 0);

oldRng = rng;
cleanupObj = onCleanup(@() rng(oldRng));
rng(seed, 'twister');
bootFH = nan(nBootstrap, 1);
for k = 1:nBootstrap
    chosen = randi(nBlocks, nBlocks, 1);
    sa = sum(blockAD(chosen));
    sb = sum(blockBC(chosen));
    if sa + sb > eps
        bootFH(k) = sa / (sa + sb);
    end
end
ci = prctile(bootFH(isfinite(bootFH)), [2.5 97.5]);

if ci(1) > 0.60 && classicFraction >= 0.60
    label = 'hyperbolic_convergence_extension_front';
    conclusion = sprintf(['The FA is classified as a hyperbolic convergence-extension front: ' ...
        'F_H=%.3f (95%% block-bootstrap CI %.3f-%.3f), and %.1f%% of valid points have a<0,d>0.'], ...
        FH, ci(1), ci(2), 100 * classicFraction);
elseif ci(2) < 0.40 && shearFraction >= 0.60
    label = 'shear_coupled_2d_filament';
    conclusion = sprintf(['The FA is classified as a shear-coupled two-dimensional filament: ' ...
        'F_H=%.3f (95%% block-bootstrap CI %.3f-%.3f), and %.1f%% of valid points have bc>0.'], ...
        FH, ci(1), ci(2), 100 * shearFraction);
else
    label = 'mixed_hyperbolic_shear_coupled';
    conclusion = sprintf(['The FA is classified as a mixed hyperbolic/shear-coupled structure: ' ...
        'F_H=%.3f (95%% block-bootstrap CI %.3f-%.3f), with %.1f%% a<0,d>0 and %.1f%% bc>0.'], ...
        FH, ci(1), ci(2), 100 * classicFraction, 100 * shearFraction);
end

classification = struct();
classification.label = label;
classification.conclusion = conclusion;
classification.F_H = FH;
classification.F_C = 1 - FH;
classification.ci95 = ci;
classification.classic_fraction = classicFraction;
classification.shear_fraction = shearFraction;
classification.n_valid = nnz(valid);
classification.n_blocks = nBlocks;
classification.n_bootstrap = nBootstrap;
classification.block_km = blockKm;
classification.bootstrap_F_H = bootFH;
end

