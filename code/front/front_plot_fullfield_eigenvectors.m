function front_plot_fullfield_eigenvectors(result, outputFile)
%FRONT_PLOT_FULLFIELD_EIGENVECTORS Plot all FL eigendirections in final style.

fig = figure('Visible', 'off', 'Color', 'w', 'Units', 'centimeters', ...
    'Position', [1 1 28 25]);
cleanupObj = onCleanup(@() close(fig));
ax = axes(fig);
pcolor(ax, result.data.lonir, result.data.latir, result.data.M2);
shading(ax, 'flat');
hold(ax, 'on');
colormap(ax, summer(256));
cb = colorbar(ax);
cb.Label.String = 'M2';

successful = result.fronts([result.fronts.success]);
hFL = plot(ax, NaN, NaN, 'k-', 'LineWidth', 1.5);
for k = 1:numel(successful)
    r = successful(k).ridge;
    plot(ax, r.x, r.y, 'k-', 'LineWidth', 1.5, 'HandleVisibility', 'off');
    mid = max(1, round(numel(r.x) / 2));
    text(ax, r.x(mid), r.y(mid), sprintf(' %d', successful(k).fa_id), ...
        'Color', 'k', 'FontSize', 7, 'FontWeight', 'bold', ...
        'BackgroundColor', 'w', 'Margin', 0.5, 'Clipping', 'on');
end

allValues = [];
for k = 1:numel(successful)
    samples = successful(k).eigen_samples;
    for j = 1:numel(samples)
        allValues = [allValues; abs(samples(j).S_values(:)); ...
            abs(samples(j).AT_values(:))]; %#ok<AGROW>
    end
end
if isempty(allValues)
    error('Front:NoEigenSamples', 'No full-field eigenvector samples are available.');
end
maxValue = max(allValues);
ratios = allValues ./ maxValue;
% Match the final ROI figure: nominal 3x maximum half-length is 6 km,
% followed by square-root compression while preserving the linear mean.
compressedMaxHalfLengthKm = 6 * mean(ratios) / mean(sqrt(ratios));
colors = {[0.68 0.68 0.68], [0.42 0.42 0.42], ...
    [0 0.25 0.95], [0.95 0.02 0.02]};

for k = 1:numel(successful)
    samples = successful(k).eigen_samples;
    for j = 1:numel(samples)
        s = samples(j);
        plot(ax, s.x, s.y, 'ko', 'MarkerFaceColor', 'w', 'MarkerSize', 2.5, ...
            'HandleVisibility', 'off');
        vectors = {s.S_vectors(:,1), s.S_vectors(:,2), ...
            s.AT_vectors(:,1), s.AT_vectors(:,2)};
        values = [s.S_values(1), s.S_values(2), s.AT_values(1), s.AT_values(2)];
        widths = [1.3 1.3 0.65 0.85];
        for q = 1:4
            len = compressedMaxHalfLengthKm * sqrt(abs(values(q)) / maxValue);
            uv = len .* vectors{q};
            draw_double_arrow(ax, s.x, s.y, uv(1), uv(2), ...
                colors{q}, widths(q), values(q) < 0);
        end
    end
end

axis(ax, 'equal');
xlim(ax, [min(result.data.lonir(:)) max(result.data.lonir(:))]);
ylim(ax, [min(result.data.latir(:)) max(result.data.latir(:))]);
set(ax, 'FontSize', 12, 'Box', 'on', 'Layer', 'top');
xlabel(ax, 'x (km)'); ylabel(ax, 'y (km)');
title(ax, sprintf(['Day 40 full-field top 10%% fronts: %d FLs, 3x arrows; ' ...
    'negative eigenvalues point inward'], numel(successful)));
hSmin = plot(ax, NaN, NaN, '-', 'Color', colors{1}, 'LineWidth', 1.3);
hSmax = plot(ax, NaN, NaN, '-', 'Color', colors{2}, 'LineWidth', 1.3);
hATneg = plot(ax, NaN, NaN, '-', 'Color', colors{3}, 'LineWidth', 0.65);
hATpos = plot(ax, NaN, NaN, '-', 'Color', colors{4}, 'LineWidth', 0.85);
legend(ax, [hFL hSmin hSmax hATneg hATpos], ...
    {'FL','S min (smaller |lambda|)','S max (larger |lambda|)', ...
    'A^T negative','A^T positive'}, 'Location', 'southoutside', ...
    'NumColumns', 3);
exportgraphics(fig, outputFile, 'Resolution', 600, 'BackgroundColor', 'white');
end

function draw_double_arrow(ax, x0, y0, ux, uy, colorValue, lineWidth, inward)
if hypot(ux, uy) <= eps, return; end
if inward
    quiver(ax, x0 + ux, y0 + uy, -ux, -uy, 0, 'Color', colorValue, ...
        'LineWidth', lineWidth, 'MaxHeadSize', 0.7, 'AutoScale', 'off', ...
        'HandleVisibility', 'off');
    quiver(ax, x0 - ux, y0 - uy, ux, uy, 0, 'Color', colorValue, ...
        'LineWidth', lineWidth, 'MaxHeadSize', 0.7, 'AutoScale', 'off', ...
        'HandleVisibility', 'off');
else
    quiver(ax, x0, y0, ux, uy, 0, 'Color', colorValue, ...
        'LineWidth', lineWidth, 'MaxHeadSize', 0.7, 'AutoScale', 'off', ...
        'HandleVisibility', 'off');
    quiver(ax, x0, y0, -ux, -uy, 0, 'Color', colorValue, ...
        'LineWidth', lineWidth, 'MaxHeadSize', 0.7, 'AutoScale', 'off', ...
        'HandleVisibility', 'off');
end
end

