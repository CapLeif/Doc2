function plot_frontogenesis_fields(data, diagnostics, outputFile, day)
%PLOT_FRONTOGENESIS_FIELDS Plot day-specific frontogenesis function and rate.

if nargin < 4, day = 40; end
mask = diagnostics.valid_mask;
functionScale = scientific_scale(diagnostics.frontogenesis_function(mask));
functionField = diagnostics.frontogenesis_function ./ functionScale;
functionLimit = robust_symmetric_limit(functionField(mask), 99);
rateLimit = robust_symmetric_limit(diagnostics.frontogenetic_rate_f(mask), 99);

fig = figure('Units', 'centimeters', 'Position', [2 2 31 13], 'Color', 'w');
t = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

ax1 = nexttile(t);
draw_field(ax1, data.lonir, data.latir, functionField, functionLimit);
title(ax1, sprintf('(a) Frontogenesis function, day %d', day), ...
    'FontWeight', 'normal');
xlabel(ax1, 'x (km)'); ylabel(ax1, 'y (km)');
cb1 = colorbar(ax1);
cb1.Label.String = sprintf('F_{hor,b} (10^{%d} s^{-5})', round(log10(functionScale)));

ax2 = nexttile(t);
draw_field(ax2, data.lonir, data.latir, diagnostics.frontogenetic_rate_f, rateLimit);
title(ax2, sprintf('(b) Frontogenetic rate, day %d', day), ...
    'FontWeight', 'normal');
xlabel(ax2, 'x (km)'); ylabel(ax2, 'y (km)');
cb2 = colorbar(ax2);
cb2.Label.String = 'T_b/f';

colormap(fig, blue_white_red(256));
exportgraphics(fig, outputFile, 'Resolution', 600, 'BackgroundColor', 'white');
close(fig);
end

function draw_field(ax, x, y, values, colorLimit)
imagesc(ax, x(1,:), y(:,1), values);
set(ax, 'YDir', 'normal', 'FontName', 'Arial', 'FontSize', 11, ...
    'Box', 'on', 'LineWidth', 1.1);
toolbar = axtoolbar(ax, {});
toolbar.Visible = 'off';
disableDefaultInteractivity(ax);
axis(ax, 'image');
caxis(ax, [-colorLimit colorLimit]);
end

function scale = scientific_scale(values)
values = abs(values(isfinite(values) & values ~= 0));
if isempty(values)
    scale = 1;
else
    scale = 10.^floor(log10(percentile(values, 99)));
end
end

function limit = robust_symmetric_limit(values, percentileValue)
values = abs(values(isfinite(values)));
if isempty(values)
    limit = 1;
else
    limit = percentile(values, percentileValue);
    if ~(isfinite(limit) && limit > 0), limit = max(values); end
    if ~(isfinite(limit) && limit > 0), limit = 1; end
end
end

function value = percentile(values, p)
values = sort(values(:));
position = 1 + (numel(values) - 1) .* p ./ 100;
lo = floor(position); hi = ceil(position);
if lo == hi
    value = values(lo);
else
    value = values(lo) + (position - lo) .* (values(hi) - values(lo));
end
end

function cmap = blue_white_red(n)
half = ceil(n / 2);
blue = [linspace(0.10,1,half).', linspace(0.25,1,half).', ones(half,1)];
red = [ones(n-half+1,1), linspace(1,0.20,n-half+1).', ...
    linspace(1,0.15,n-half+1).'];
cmap = [blue; red(2:end,:)];
end

