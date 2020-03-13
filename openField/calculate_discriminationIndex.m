function DI = calculate_discriminationIndex(sumFiringRateObjectAll)
DI = nan(size(sumFiringRateObjectAll{1},1),length(sumFiringRateObjectAll));
for i = 1:length(sumFiringRateObjectAll)
    if size(sumFiringRateObjectAll{i},1) == 1
        if ~isnan(sumFiringRateObjectAll{i})
            DI(1,i) = diff(sumFiringRateObjectAll{i})/sum(sumFiringRateObjectAll{i});
        end
    else
        if sum(sumFiringRateObjectAll{i}(:)) ~= 0
            DI(:,i) = (sumFiringRateObjectAll{i}(:,2)-sumFiringRateObjectAll{i}(:,1))./sum(sumFiringRateObjectAll{i},2);
        end
    end
end
