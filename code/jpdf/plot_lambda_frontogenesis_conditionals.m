function plot_lambda_frontogenesis_conditionals(stats, outputFile, day, statistic)
%PLOT_LAMBDA_FRONTOGENESIS_CONDITIONALS Plot conditional mean or median on lambda-JPDF.

if nargin < 3, day = 40; end
if nargin < 4 || isempty(statistic), statistic = 'mean'; end
switch lower(statistic)
    case 'mean'
        frontogenesis = stats.conditional_frontogenesis_function;
        rate_f = stats.conditional_frontogenetic_rate_f;
        operator = 'E';
        titleWord = 'mean';
    case 'median'
        frontogenesis = stats.conditional_median_frontogenesis_function;
        rate_f = stats.conditional_median_frontogenetic_rate_f;
        operator = 'median';
        titleWord = 'median';
    otherwise
        error('Frontogenesis:InvalidConditionalStatistic', ...
            'statistic must be either mean or median.');
end
functionScale = scientific_scale(frontogenesis(:));
frontogenesis = frontogenesis ./ functionScale;

fig = figure('Units', 'centimeters', 'Position', [2 2 32 14], 'Color', 'w');
t = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

ax1 = nexttile(t);
draw_conditional(ax1, stats, frontogenesis, ...
    sprintf('%s[F_{hor,b}|lambda-JPDF] (10^{%d} s^{-5})', operator, ...
    round(log10(functionScale))));
title(ax1, sprintf('(a) Conditional %s frontogenesis function, day %d', ...
    titleWord, day), ...
    'FontWeight', 'normal');

ax2 = nexttile(t);
draw_conditional(ax2, stats, rate_f, ...
    sprintf('%s[T_b/f|lambda-JPDF]', operator));
title(ax2, sprintf('(b) Conditional %s frontogenetic rate, day %d', ...
    titleWord, day), ...
    'FontWeight', 'normal');

colormap(fig, blue_white_red(256));
sgtitle(t, sprintf('|Ro| <= 10; bin = %.2f; minimum %d samples/bin; N = %d', ...
    diff(stats.bin_edges(1:2)), stats.min_count, stats.in_range_count), ...
    'FontWeight', 'normal', 'FontSize', 11);
exportgraphics(fig, outputFile, 'Resolution', 600, 'BackgroundColor', 'white');
close(fig);
end

function draw_conditional(ax, stats, values, colorbarLabel)
centers = stats.bin_centers;
hm = imagesc(ax, centers, centers, values.');
set(hm, 'AlphaData', isfinite(values.'));
set(ax, 'YDir', 'normal', 'FontName', 'Arial', 'FontSize', 11, ...
    'Box', 'on', 'LineWidth', 1.1, 'Color', [0.93 0.93 0.93]);
toolbar = axtoolbar(ax, {});
toolbar.Visible = 'off';
disableDefaultInteractivity(ax);
hold(ax, 'on');
limit = robust_symmetric_limit(values(:), 98);
caxis(ax, [-limit limit]);
cb = colorbar(ax);
cb.Label.String = colorbarLabel;

[X, Y] = meshgrid(centers, centers);
logPdf = log10(stats.pdf.');
pdfLevels = [-3 -2 -1];
contour(ax, X, Y, logPdf, pdfLevels, 'Color', [0.15 0.15 0.15], ...
    'LineWidth', 0.7, 'LineStyle', ':');

lambdaMinus = 0.5 .* (X + Y);
lambdaMinus(Y >= 0) = NaN;
[C, h] = contour(ax, X, Y, lambdaMinus, [-4 -3 -2 -1 -0.5], ...
    'Color', [0.45 0.45 0.45], 'LineWidth', 1.0, 'LineStyle', '--');
clabel(C, h, 'FontSize', 8, 'Color', [0.35 0.35 0.35], ...
    'LabelSpacing', 300);

lims = [stats.bin_edges(1) stats.bin_edges(end)];
plot(ax, lims, lims, 'k--', 'LineWidth', 0.8);
plot(ax, lims, -lims, 'k--', 'LineWidth', 0.8);
xline(ax, 0, 'k--', 'LineWidth', 0.8);
yline(ax, 0, 'k--', 'LineWidth', 0.8);
xlim(ax, lims); ylim(ax, lims); axis(ax, 'square');
xlabel(ax, '$\Delta/f$', 'Interpreter', 'latex');
ylabel(ax, '$-\sqrt{\Omega}/f\ (\Omega\geq0);\ +\sqrt{-\Omega}/f\ (\Omega<0)$', ...
    'Interpreter', 'latex');
text(ax, 0.03, 0.03, 'dashed gray: \lambda_-/f; dotted black: log_{10} PDF', ...
    'Units', 'normalized', 'FontSize', 8, 'Color', [0.2 0.2 0.2], ...
    'BackgroundColor', 'w', 'Margin', 2);
end

function scale = scientific_scale(values)
values = abs(values(isfinite(values) & values ~= 0));
if isempty(values)
    scale = 1;
else
    scale = 10.^floor(log10(percentile(values, 98)));
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

