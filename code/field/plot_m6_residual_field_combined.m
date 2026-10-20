function plot_m6_residual_field_combined(result, outputFile, fullMask, filamentMask, colorLim)
%PLOT_M6_RESIDUAL_FIELD_COMBINED Plot full-field and filament det residual maps.

if nargin < 3 || isempty(fullMask)
    fullMask = true(size(result.residual_f2));
end
if nargin < 4 || isempty(filamentMask)
    filamentMask = result.filament_mask;
end
if nargin < 5 || isempty(colorLim)
    colorLim = [-1 1];
end
if ~isequal(size(fullMask), size(result.residual_f2)) || ...
        ~isequal(size(filamentMask), size(result.residual_f2))
    error('M6:InvalidFieldMask', 'Masks must match the residual field size.');
end

residual_map = result.residual_f2 ./ (result.Strain_f.^2);
finiteMask = isfinite(residual_map) & result.Strain_f.^2 > eps;
fullMask = fullMask & finiteMask;
filamentMask = filamentMask & finiteMask;

if nnz(fullMask) == 0 || nnz(filamentMask) == 0
    error('M6:EmptyResidualField', 'No finite residual values remain after masking.');
end

edges = colorLim(1):0.1:colorLim(2);
if abs(edges(end) - colorLim(2)) > 10 * eps(colorLim(2))
    error('M6:InvalidColorLim', 'colorLim must align with 0.1 intervals.');
end
cmap = diverging_colormap(numel(edges) - 1);

fig = figure('Units', 'centimeters', 'Position', [2 2 42 34], 'Color', 'w');

ax1 = axes(fig, 'Position', [0.040 0.515 0.450 0.405]);
plot_panel(ax1, result, residual_map, fullMask, cmap, colorLim, ...
    sprintf('(a) Full field  Day %d  N=%d', result.day, nnz(fullMask)));

ax2 = axes(fig, 'Position', [0.510 0.515 0.450 0.405]);
plot_panel(ax2, result, residual_map, filamentMask, cmap, colorLim, ...
    sprintf('(b) Filaments  Day %d  N=%d', result.day, nnz(filamentMask)));

draw_discrete_colorbar(fig, cmap, edges, '$-4\det(\Lambda)/\sigma^2$');

ax3 = axes(fig, 'Position', [0.105 0.070 0.835 0.245]);
plot_theta135_pdf_panel(ax3, result, residual_map, filamentMask, 0.2);
drawnow;
exportgraphics(fig, outputFile, 'Resolution', 600, 'BackgroundColor', 'white');
close(fig);
end

function plot_panel(ax, result, residual_map, sampleMask, cmap, colorLim, titleText)
panel_map = residual_map;
panel_map(~sampleMask) = NaN;

pcolor(ax, result.lonir, result.latir, panel_map);
shading(ax, 'flat');
axis(ax, 'equal');
xlim(ax, [0 200]);
ylim(ax, [0 200]);
colormap(ax, cmap);
caxis(ax, colorLim);
set(ax, 'Color', [0.82 0.82 0.82], 'FontSize', 18, 'Box', 'on', ...
    'LineWidth', 1.6, 'FontName', 'Arial', 'TickDir', 'out', ...
    'XTick', 0:50:200, 'YTick', 0:50:200);
xlabel(ax, '$x$', 'Interpreter', 'latex', 'FontSize', 21);
ylabel(ax, '$y$', 'Interpreter', 'latex', 'FontSize', 21);
titleHandle = title(ax, titleText, 'FontWeight', 'normal', 'FontSize', 22, ...
    'Interpreter', 'none');
set(titleHandle, 'Units', 'normalized', 'Position', [0.5 1.045 0]);
end

function draw_discrete_colorbar(fig, cmap, edges, labelText)
cbAx = axes(fig, 'Position', [0.12 0.425 0.76 0.024]);
hold(cbAx, 'on');
for k = 1:size(cmap, 1)
    patch(cbAx, [edges(k) edges(k + 1) edges(k + 1) edges(k)], ...
        [0 0 1 1], cmap(k, :), 'EdgeColor', [0.20 0.20 0.20], ...
        'LineWidth', 0.45);
end
plot(cbAx, [0 0], [0 1], 'k-', 'LineWidth', 1.5);
xlim(cbAx, [edges(1) edges(end)]);
ylim(cbAx, [0 1]);
labelTicks = edges(1):0.2:edges(end);
set(cbAx, 'YTick', [], 'XTick', labelTicks, 'FontSize', 12, ...
    'FontName', 'Arial', 'Box', 'on', 'LineWidth', 1.2, ...
    'TickDir', 'out', 'TickLength', [0.004 0.004]);
cbAx.XTickLabel = compose('%.1f', labelTicks);
cbAx.XTickLabelRotation = 0;
title(cbAx, labelText, 'Interpreter', 'latex', 'FontSize', 16, ...
    'FontWeight', 'normal');
end

function plot_theta135_pdf_panel(ax, result, residual, sampleMask, binWidth)
baseValid = isfinite(result.theta_all) & isfinite(residual) & ...
    result.Strain_f.^2 > eps & sampleMask;
edges = -1:binWidth:1;
centers = (edges(1:end-1) + edges(2:end)) ./ 2;
halfWidths = [10 5 1];
colors = [0.00 0.45 0.74; 0.85 0.33 0.10; 0.47 0.67 0.19];
hold(ax, 'on');
legendHandles = gobjects(numel(halfWidths), 1);
legendLabels = cell(numel(halfWidths), 1);
for k = 1:numel(halfWidths)
    width = halfWidths(k);
    mask = baseValid & result.theta_all >= 135 - width & ...
        result.theta_all <= 135 + width & residual >= edges(1) & ...
        residual <= edges(end);
    values = residual(mask);
    if isempty(values)
        pdfValues = zeros(size(centers));
    else
        pdfValues = histcounts(values, edges, 'Normalization', 'pdf');
    end
    legendHandles(k) = plot(ax, centers, pdfValues, 'Color', colors(k, :), ...
        'LineStyle', '-', 'LineWidth', 3.0);
    legendLabels{k} = sprintf('$\\pm%d^\\circ$ ($N=%d$)', width, numel(values));
end
xline(ax, 0, 'Color', [0.55 0.55 0.55], 'LineStyle', ':', ...
    'LineWidth', 1.8, 'HandleVisibility', 'off');
grid(ax, 'on');
box(ax, 'on');
xlim(ax, [-1 1]);
currentYLim = ylim(ax);
ylim(ax, [0 currentYLim(2) * 1.08]);
set(ax, 'FontName', 'Arial', 'FontSize', 17, 'LineWidth', 1.5, ...
    'TickDir', 'out', 'XTick', -1:0.2:1);
xlabel(ax, '$-4\det(\Lambda)/\sigma^2$', 'Interpreter', 'latex', ...
    'FontSize', 21);
ylabel(ax, 'PDF', 'FontSize', 20);
titleHandle = title(ax, '(c) Residual PDF near \theta=135^\circ  filaments  bin=0.2', ...
    'FontWeight', 'normal', 'FontSize', 22, 'Interpreter', 'tex');
set(titleHandle, 'Units', 'normalized', 'Position', [0.5 1.035 0]);
legend(ax, legendHandles, legendLabels, 'Location', 'northeast', ...
    'FontSize', 16, 'Box', 'off', 'Interpreter', 'latex');
end

function cmap = diverging_colormap(n)
blue = [0.10 0.25 0.65];
white = [1.00 1.00 1.00];
red = [0.70 0.10 0.10];
half = ceil(n / 2);
c1 = [linspace(blue(1), white(1), half)', ...
      linspace(blue(2), white(2), half)', ...
      linspace(blue(3), white(3), half)'];
c2 = [linspace(white(1), red(1), n - half)', ...
      linspace(white(2), red(2), n - half)', ...
      linspace(white(3), red(3), n - half)'];
cmap = [c1; c2];
end

