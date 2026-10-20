function samples = front_sample_eigenvectors(ridge, data, searchKm, fractions)
%FRONT_SAMPLE_EIGENVECTORS Sample S and A' eigensystems along FL.

if nargin < 3 || isempty(searchKm)
    searchKm = 1;
end
if nargin < 4 || isempty(fractions)
    fractions = [0 0.25 0.5 0.75 1];
end
fractions = fractions(:).';
if any(~isfinite(fractions)) || any(fractions < 0 | fractions > 1) || ...
        any(diff(fractions) <= 0)
    error('Front:InvalidFractions', ...
        'fractions must be finite, strictly increasing values from 0 to 1.');
end
targetArc = fractions .* ridge.length_km;
nSample = numel(fractions);
samples = repmat(struct('fraction', [], 'target_arc_km', [], 'arc_km', [], ...
    'x', [], 'y', [], 'A', [], 'S', [], 'S_values', [], 'S_vectors', [], ...
    'AT_values', [], 'AT_vectors', [], 'relocated_km', []), 1, nSample);

for k = 1:nSample
    distanceToTarget = abs(ridge.arc_km - targetArc(k));
    candidates = find(distanceToTarget <= searchKm + eps);
    [~, order] = sort(distanceToTarget(candidates));
    candidates = candidates(order);
    found = false;
    for q = 1:numel(candidates)
        j = candidates(q);
        A = interpolate_A(data, ridge.x(j), ridge.y(j));
        [VAT, DAT] = eig(A.');
        evalAT = diag(DAT);
        if max(abs(imag(evalAT))) <= 1e-10 && min(real(evalAT)) < 0 && max(real(evalAT)) > 0
            evalAT = real(evalAT);
            VAT = real(VAT);
            [evalAT, indAT] = sort(evalAT, 'ascend');
            VAT = VAT(:, indAT);
            S = 0.5 .* (A + A.');
            [VS, DS] = eig(S);
            evalS = diag(DS);
            % S min/max denote the smaller/larger eigenvalue magnitude,
            % not the algebraically smaller/larger signed eigenvalue.
            [~, indS] = sort(abs(evalS), 'ascend');
            evalS = evalS(indS);
            VS = VS(:, indS);
            samples(k).fraction = fractions(k);
            samples(k).target_arc_km = targetArc(k);
            samples(k).arc_km = ridge.arc_km(j);
            samples(k).x = ridge.x(j);
            samples(k).y = ridge.y(j);
            samples(k).A = A;
            samples(k).S = S;
            samples(k).S_values = evalS;
            samples(k).S_vectors = VS;
            samples(k).AT_values = evalAT;
            samples(k).AT_vectors = VAT;
            samples(k).relocated_km = distanceToTarget(j);
            found = true;
            break;
        end
    end
    if ~found
        error('Front:NoHyperbolicSample', ...
            'No real negative-positive A'' eigenpair exists within %.1f km of FL fraction %.2f.', ...
            searchKm, fractions(k));
    end
end

for k = 2:nSample
    for mode = 1:2
        if dot(samples(k - 1).S_vectors(:, mode), samples(k).S_vectors(:, mode)) < 0
            samples(k).S_vectors(:, mode) = -samples(k).S_vectors(:, mode);
        end
        if dot(samples(k - 1).AT_vectors(:, mode), samples(k).AT_vectors(:, mode)) < 0
            samples(k).AT_vectors(:, mode) = -samples(k).AT_vectors(:, mode);
        end
    end
end
end

function A = interpolate_A(data, xq, yq)
ux = interp2(data.lonir, data.latir, data.pupx, xq, yq, 'linear');
uy = interp2(data.lonir, data.latir, data.pupy, xq, yq, 'linear');
vx = interp2(data.lonir, data.latir, data.pvpx, xq, yq, 'linear');
vy = interp2(data.lonir, data.latir, data.pvpy, xq, yq, 'linear');
A = [ux uy; vx vy];
if any(~isfinite(A(:)))
    error('Front:InterpolationFailed', 'Velocity-gradient interpolation returned nonfinite values.');
end
end

