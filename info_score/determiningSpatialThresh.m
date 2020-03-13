function threshSpatial = determiningSpatialThresh(FiringRateAll,segments)
if ~exist('segment','var') || isempty(segments)
    segments = 1:size(FiringRateAll,1);
end
k = 0;
maxCount = zeros(length(segments)*size(FiringRateAll,2),1);
for sessionIndex = 1:size(FiringRateAll,2)
    firingRate = FiringRateAll(:,sessionIndex);
    for i = segments
        k = k+1;
        if ~isempty(firingRate{i})
            maxCount(k) = max(firingRate{i}(:));
        end
    end
end

threshSpatial = quantile(maxCount,0.85);