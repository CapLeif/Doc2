function discovery = front_extract_fullfield_fas(data, percentage, minMajorAxisKm)
%FRONT_EXTRACT_FULLFIELD_FAS Discover elongated top-M2 connected components.

if nargin < 2 || isempty(percentage), percentage = 10; end
if nargin < 3 || isempty(minMajorAxisKm), minMajorAxisKm = 10; end
required = {'M2','lonir','latir'};
for k = 1:numel(required)
    if ~isfield(data, required{k})
        error('Front:MissingVariable', 'Required field missing: %s', required{k});
    end
end
if ~(isscalar(percentage) && percentage > 0 && percentage < 100)
    error('Front:InvalidPercentage', 'percentage must lie strictly between 0 and 100.');
end

valid = isfinite(data.M2) & isfinite(data.lonir) & isfinite(data.latir);
threshold = prctile(data.M2(valid), 100 - percentage);
rawMask = valid & data.M2 >= threshold;
cc = bwconncomp(rawMask, 8);
props = regionprops(cc, 'Area', 'MajorAxisLength', 'Eccentricity', ...
    'Centroid', 'BoundingBox', 'PixelIdxList');
dx = median(abs(diff(data.lonir, 1, 2)), 'all', 'omitnan');
dy = median(abs(diff(data.latir, 1, 1)), 'all', 'omitnan');
scaleKm = sqrt(dx * dy);

n = cc.NumObjects;
componentId = (1:n).';
pixelCount = reshape([props.Area], [], 1);
areaKm2 = pixelCount .* dx .* dy;
majorAxisKm = reshape([props.MajorAxisLength], [], 1) .* scaleKm;
eccentricity = reshape([props.Eccentricity], [], 1);
centroids = reshape([props.Centroid], 2, []).';
centroidXKm = interp2(data.lonir, centroids(:,1), centroids(:,2), 'linear');
centroidYKm = interp2(data.latir, centroids(:,1), centroids(:,2), 'linear');
accepted = majorAxisKm >= minMajorAxisKm;

acceptedIds = componentId(accepted);
[~, order] = sortrows([-centroidYKm(accepted), centroidXKm(accepted)], [1 2]);
acceptedIds = acceptedIds(order);

discovery = struct();
discovery.percentage = percentage;
discovery.threshold = threshold;
discovery.raw_mask = rawMask;
discovery.component_count_raw = n;
discovery.component_count_accepted = numel(acceptedIds);
discovery.min_major_axis_km = minMajorAxisKm;
discovery.dx_km = dx;
discovery.dy_km = dy;
discovery.cc = cc;
discovery.props = props;
discovery.accepted_component_ids = acceptedIds;
discovery.table = table(componentId, pixelCount, areaKm2, majorAxisKm, ...
    eccentricity, centroidXKm, centroidYKm, accepted, ...
    'VariableNames', {'component_id','pixel_count','area_km2', ...
    'major_axis_km','eccentricity','centroid_x_km','centroid_y_km','accepted_FA'});
end

