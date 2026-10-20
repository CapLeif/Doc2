function [frontsOut, audit] = front_plan_fullfield_eigen_samples(fronts, data, referenceFlId)
%FRONT_PLAN_FULLFIELD_EIGEN_SAMPLES Assign length-aware arrow samples.

if nargin < 3 || isempty(referenceFlId), referenceFlId = 44; end
referenceIndex = find([fronts.fl_id] == referenceFlId, 1);
if isempty(referenceIndex)
    error('Front:MissingReferenceFL', 'Reference FL%d was not found.', referenceFlId);
end
referenceLength = fronts(referenceIndex).ridge.length_km;
frontsOut = fronts;
audit = table();

for k = 1:numel(fronts)
    ridge = fronts(k).ridge;
    L = ridge.length_km;
    if fronts(k).fl_id == referenceFlId
        targetArc = [0 L/2 L];
        rule = "reference_end_center_end";
    elseif L < referenceLength
        targetArc = L/2;
        rule = "shorter_than_FL44_center_only";
    else
        nIntervals = max(2, round(L/10));
        targetArc = linspace(0, L, nIntervals+1);
        rule = "approximately_10km_including_endpoints";
    end
    samples = [];
    failures = strings(0,1);
    completeCount = 0;
    for q = 1:numel(targetArc)
        fraction = targetArc(q) / L;
        try
            s = front_sample_eigenvectors(ridge, data, 1, fraction);
            s.planned_arc_km = targetArc(q);
            samples = [samples s]; %#ok<AGROW>
            completeCount = completeCount + 1;
        catch ME
            if strcmp(ME.identifier, 'Front:NoHyperbolicSample')
                failures(end+1,1) = sprintf('arc %.3f km: %s', targetArc(q), ME.message); %#ok<AGROW>
                s = sample_s_only(ridge, data, targetArc(q), fraction);
                samples = [samples s]; %#ok<AGROW>
            else
                rethrow(ME);
            end
        end
    end
    frontsOut(k).eigen_samples = samples;
    frontsOut(k).arrow_target_arc_km = targetArc;
    frontsOut(k).arrow_rule = rule;
    frontsOut(k).arrow_failures = failures;
    audit = [audit; table(fronts(k).fl_id, L, referenceLength, rule, ...
        numel(targetArc), completeCount, numel(targetArc)-completeCount, ...
        strjoin(failures, ' | '), ...
        'VariableNames', {'fl_id','FL_length_km','FL44_length_km','sampling_rule', ...
        'planned_arrow_points','complete_four_arrow_points','S_only_points', ...
        'failure_details'})]; %#ok<AGROW>
end
end

function s = sample_s_only(ridge, data, targetArc, fraction)
xq = interp1(ridge.arc_km, ridge.x, targetArc, 'linear');
yq = interp1(ridge.arc_km, ridge.y, targetArc, 'linear');
ux = interp2(data.lonir,data.latir,data.pupx,xq,yq,'linear');
uy = interp2(data.lonir,data.latir,data.pupy,xq,yq,'linear');
vx = interp2(data.lonir,data.latir,data.pvpx,xq,yq,'linear');
vy = interp2(data.lonir,data.latir,data.pvpy,xq,yq,'linear');
A = [ux uy; vx vy];
S = 0.5.*(A+A.');
[VS,DS] = eig(S);
evalS = diag(DS);
[~,order] = sort(abs(evalS),'ascend');
s = struct('fraction',fraction,'target_arc_km',targetArc,'arc_km',targetArc, ...
    'x',xq,'y',yq,'A',A,'S',S,'S_values',evalS(order), ...
    'S_vectors',VS(:,order),'AT_values',[NaN;NaN], ...
    'AT_vectors',nan(2,2),'relocated_km',0,'planned_arc_km',targetArc);
end

