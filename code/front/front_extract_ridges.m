function ridges = front_extract_ridges(M2, pbhpx, pbhpy, x, y, faMask, sigmaGrid, minLengthKm, maxTangentAngleDeg)
%FRONT_EXTRACT_RIDGES Extract multiple skeleton-separated constrained ridges.

if nargin < 7 || isempty(sigmaGrid), sigmaGrid = 1; end
if nargin < 8 || isempty(minLengthKm), minLengthKm = 8; end
if nargin < 9 || isempty(maxTangentAngleDeg), maxTangentAngleDeg = 45; end
if ~isequal(size(M2), size(pbhpx), size(pbhpy), size(x), size(y), size(faMask))
    error('Front:SizeMismatch', 'Ridge inputs must have equal sizes.');
end

dx = median(abs(diff(x, 1, 2)), 'all', 'omitnan');
dy = median(abs(diff(y, 1, 1)), 'all', 'omitnan');
h = min(dx, dy);
M2filled = fillmissing(M2, 'nearest', 1);
M2filled = fillmissing(M2filled, 'nearest', 2);
M2smooth = imgaussfilt(M2filled, sigmaGrid, 'FilterSize', 5, 'Padding', 'replicate');
[gx, gy] = gradient(M2smooth, dx, dy);
[gxx, gxy1] = gradient(gx, dx, dy);
[gxy2, gyy] = gradient(gy, dx, dy);
gxy = 0.5 .* (gxy1 + gxy2);

gradB = hypot(pbhpx, pbhpy);
nx = pbhpx ./ gradB; ny = pbhpy ./ gradB;
tx = -ny; ty = nx;
normalDerivative = nx .* gx + ny .* gy;
normalCurvature = nx.^2 .* gxx + 2 .* nx .* ny .* gxy + ny.^2 .* gyy;
offset = -normalDerivative ./ normalCurvature;
candidateMask = faMask & isfinite(offset) & isfinite(normalCurvature) & ...
    normalCurvature < 0 & abs(offset) <= 0.5 * hypot(dx, dy) & ...
    isfinite(nx) & isfinite(ny);

candidateIndex = find(candidateMask);
xcGrid = nan(size(M2)); ycGrid = nan(size(M2));
xcGrid(candidateIndex) = x(candidateIndex) + offset(candidateIndex).*nx(candidateIndex);
ycGrid(candidateIndex) = y(candidateIndex) + offset(candidateIndex).*ny(candidateIndex);
inside = interp2(x, y, double(faMask), xcGrid(candidateIndex), ...
    ycGrid(candidateIndex), 'nearest', 0) > 0.5;
validMask = false(size(candidateMask));
validMask(candidateIndex(inside)) = true;
if nnz(validMask) < 5
    ridges = empty_ridges(); return;
end

% Thin the strict physical ridge support before graph construction. This
% prevents a broad candidate band from creating shortcuts across separate
% fronts. Remove one grid cell around junctions so branches become distinct.
thinMask = bwmorph(validMask, 'thin', Inf);
branchMask = bwmorph(thinMask, 'branchpoints');
if any(branchMask(:))
    splitMask = thinMask & ~imdilate(branchMask, strel('square', 3));
else
    splitMask = thinMask;
end
cc = bwconncomp(splitMask, 8);
ridges = empty_ridges();
for component = 1:cc.NumObjects
    pixels = cc.PixelIdxList{component};
    if numel(pixels) < 5, continue; end
    [rows, cols] = ind2sub(size(splitMask), pixels);
    G = pixel_graph(rows, cols, xcGrid(pixels), ycGrid(pixels), ...
        tx(pixels), ty(pixels), size(splitMask), maxTangentAngleDeg);
    if numedges(G) == 0, continue; end
    bins = conncomp(G);
    for group = unique(bins)
        nodes = find(bins == group);
        if numel(nodes) < 5, continue; end
        subG = subgraph(G, nodes);
        [localPath, graphLength] = approximate_diameter(subG);
        if graphLength < minLengthKm, continue; end
        pathPixels = pixels(nodes(localPath));
        pathX = xcGrid(pathPixels); pathY = ycGrid(pathPixels);
        [ridgeX, ridgeY, sampleArc, used] = resample_path(pathX, pathY, h);
        if sampleArc(end) < minLengthKm, continue; end
        insidePath = interp2(x, y, double(faMask), ridgeX, ridgeY, 'nearest', 0) > 0.5;
        curvaturePath = interp2(x, y, normalCurvature, ridgeX, ridgeY, 'linear');
        validPoint = insidePath & curvaturePath < 0;
        % Do not discard a whole long ridge because one interpolated point
        % fails. Split at invalid points and retain every valid >=8 km run.
        runs = contiguous_true_runs(validPoint);
        for runIndex = 1:size(runs,1)
            q1 = runs(runIndex,1); q2 = runs(runIndex,2);
            segmentX = ridgeX(q1:q2); segmentY = ridgeY(q1:q2);
            segmentArc = [0; cumsum(hypot(diff(segmentX),diff(segmentY)))];
            if segmentArc(end) < minLengthKm, continue; end
            ridge = struct('x',segmentX,'y',segmentY,'arc_km',segmentArc, ...
                'length_km',segmentArc(end),'sigma_grid',sigmaGrid, ...
                'raw_candidate_indices',pathPixels(used),'dx_km',dx,'dy_km',dy, ...
                'fraction_inside_FA',1,'fraction_negative_normal_curvature',1);
            ridges(end+1,1) = ridge; %#ok<AGROW>
        end
    end
end
if isempty(ridges), return; end
[~, order] = sort([ridges.length_km], 'descend');
ridges = ridges(order);
end

function G = pixel_graph(rows, cols, xc, yc, tx, ty, maskSize, maxTangentAngleDeg)
n = numel(rows);
nodeMap = zeros(maskSize);
nodeMap(sub2ind(maskSize, rows, cols)) = 1:n;
edgeI = []; edgeJ = []; edgeW = [];
offsets = [0 1; 1 -1; 1 0; 1 1];
minAlignment = cosd(maxTangentAngleDeg);
for q = 1:size(offsets,1)
    rr = rows + offsets(q,1); cc = cols + offsets(q,2);
    valid = rr>=1 & rr<=maskSize(1) & cc>=1 & cc<=maskSize(2);
    source = find(valid);
    target = nodeMap(sub2ind(maskSize, rr(valid), cc(valid)));
    keep = target > 0;
    source = source(keep); target = target(keep);
    if isempty(source), continue; end
    dx = xc(target) - xc(source);
    dy = yc(target) - yc(source);
    stepLength = hypot(dx, dy);
    txSource = tx(source); tySource = ty(source);
    txTarget = tx(target); tyTarget = ty(target);
    normSource = hypot(txSource, tySource);
    normTarget = hypot(txTarget, tyTarget);
    alignmentSource = abs((dx .* txSource + dy .* tySource) ./ ...
        max(eps, stepLength .* normSource));
    alignmentTarget = abs((dx .* txTarget + dy .* tyTarget) ./ ...
        max(eps, stepLength .* normTarget));
    keepEdge = isfinite(stepLength) & stepLength > 0 & ...
        isfinite(alignmentSource) & isfinite(alignmentTarget) & ...
        alignmentSource >= minAlignment & alignmentTarget >= minAlignment;
    source = source(keepEdge); target = target(keepEdge); stepLength = stepLength(keepEdge);
    edgeI = [edgeI; source]; %#ok<AGROW>
    edgeJ = [edgeJ; target]; %#ok<AGROW>
    edgeW = [edgeW; stepLength]; %#ok<AGROW>
end
if isempty(edgeI), G=graph([],[],[],n); else, G=graph(edgeI,edgeJ,edgeW,n); end
end

function [path, lengthValue] = approximate_diameter(G)
d0 = distances(G,1); [~,a] = max(d0);
d1 = distances(G,a); [lengthValue,b] = max(d1);
path = shortestpath(G,a,b);
end

function [ridgeX, ridgeY, sampleArc, used] = resample_path(pathX, pathY, h)
if pathY(1)>pathY(end) || (pathY(1)==pathY(end) && pathX(1)>pathX(end))
    pathX=flip(pathX); pathY=flip(pathY);
end
arc=[0;cumsum(hypot(diff(pathX),diff(pathY)))];
[arc,used]=unique(arc,'stable'); pathX=pathX(used); pathY=pathY(used);
sampleArc=(0:h:arc(end)).';
if sampleArc(end)<arc(end), sampleArc(end+1,1)=arc(end); end
ridgeX=interp1(arc,pathX,sampleArc,'linear');
ridgeY=interp1(arc,pathY,sampleArc,'linear');
end

function runs = contiguous_true_runs(mask)
mask = mask(:);
edges = diff([false; mask; false]);
starts = find(edges == 1);
stops = find(edges == -1) - 1;
runs = [starts stops];
end

function ridges = empty_ridges()
ridges=struct('x',{},'y',{},'arc_km',{},'length_km',{},'sigma_grid',{}, ...
    'raw_candidate_indices',{},'dx_km',{},'dy_km',{}, ...
    'fraction_inside_FA',{},'fraction_negative_normal_curvature',{});
end

