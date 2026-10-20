function masks = field_filter_mask(baseMask, Ro, lambdaMinusOverF)
%FIELD_FILTER_MASK Build the shared finite, |Ro|, and lambda-minus masks.
%   BASEMASK must be supplied by the caller from the authoritative spatial
%   domain and flow-partition definition. This function does not infer VGL
%   identifiers or construct a non-SD mask from numeric category codes.

if ~isequal(size(baseMask), size(Ro), size(lambdaMinusOverF))
    error('Doc2:MaskSizeMismatch', 'All masks and fields must have equal sizes.');
end
baseMask = logical(baseMask) & isfinite(Ro) & isfinite(lambdaMinusOverF);
roMask = abs(Ro) >= 1;
lambdaMask = lambdaMinusOverF < -1;

masks = struct();
masks.base = baseMask;
masks.ro = baseMask & roMask;
masks.lambda = baseMask & lambdaMask;
masks.ro_lambda = baseMask & roMask & lambdaMask;
end

