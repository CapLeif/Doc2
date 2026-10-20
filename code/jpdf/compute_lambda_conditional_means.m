function stats = compute_lambda_conditional_means(diagnostics, sampleMask, binEdges, minCount)
%COMPUTE_LAMBDA_CONDITIONAL_MEANS Bin conditional means and medians on lambda-JPDF axes.

if nargin < 2 || isempty(sampleMask)
    sampleMask = true(size(diagnostics.delta_f));
end
if nargin < 3 || isempty(binEdges)
    binEdges = -5:0.05:5;
end
if nargin < 4 || isempty(minCount)
    minCount = 20;
end
if ~isequal(size(sampleMask), size(diagnostics.delta_f))
    error('Frontogenesis:MaskSizeMismatch', ...
        'sampleMask must match the diagnostic field size.');
end
if numel(binEdges) < 3 || any(~isfinite(binEdges)) || any(diff(binEdges) <= 0)
    error('Frontogenesis:InvalidBinEdges', ...
        'binEdges must be a finite, strictly increasing vector.');
end
if ~(isscalar(minCount) && minCount >= 1 && minCount == floor(minCount))
    error('Frontogenesis:InvalidMinCount', ...
        'minCount must be a positive integer.');
end

valid = sampleMask & diagnostics.valid_mask;
x = diagnostics.delta_f(valid);
y = diagnostics.lambda_jpdf_y(valid);
frontogenesis = diagnostics.frontogenesis_function(valid);
rate_f = diagnostics.frontogenetic_rate_f(valid);
x = x(:);
y = y(:);
frontogenesis = frontogenesis(:);
rate_f = rate_f(:);

[counts, ~, ~, xBin, yBin] = histcounts2(x, y, binEdges, binEdges);
inside = xBin > 0 & yBin > 0 & isfinite(frontogenesis) & isfinite(rate_f);
nBins = numel(binEdges) - 1;
subscripts = [xBin(inside), yBin(inside)];
sumFrontogenesis = accumarray(subscripts, frontogenesis(inside), ...
    [nBins nBins], @sum, 0);
sumRate = accumarray(subscripts, rate_f(inside), [nBins nBins], @sum, 0);

conditionalFrontogenesis = sumFrontogenesis ./ counts;
conditionalRate = sumRate ./ counts;
medianFrontogenesis = accumarray(subscripts, frontogenesis(inside), ...
    [nBins nBins], @median, NaN);
medianRate = accumarray(subscripts, rate_f(inside), ...
    [nBins nBins], @median, NaN);
enough = counts >= minCount;
conditionalFrontogenesis(~enough) = NaN;
conditionalRate(~enough) = NaN;
medianFrontogenesis(~enough) = NaN;
medianRate(~enough) = NaN;

binWidth = diff(binEdges);
if max(abs(binWidth - binWidth(1))) > 10 * eps(max(abs(binEdges)))
    error('Frontogenesis:NonuniformBins', ...
        'The lambda-JPDF implementation requires uniform bins.');
end
pdfValues = counts ./ (sum(counts(:)) .* binWidth(1).^2);
pdfValues(counts == 0) = NaN;

stats = struct();
stats.bin_edges = binEdges;
stats.bin_centers = 0.5 .* (binEdges(1:end-1) + binEdges(2:end));
stats.counts = counts;
stats.pdf = pdfValues;
stats.conditional_frontogenesis_function = conditionalFrontogenesis;
stats.conditional_frontogenetic_rate_f = conditionalRate;
stats.conditional_median_frontogenesis_function = medianFrontogenesis;
stats.conditional_median_frontogenetic_rate_f = medianRate;
stats.min_count = minCount;
stats.sample_count = nnz(valid);
stats.in_range_count = sum(counts(:));
end

