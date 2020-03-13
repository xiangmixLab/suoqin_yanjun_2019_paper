function [infoPerSecond, infoPerSpike] = comparisonSpatialInfo(firingMap, MeanFiringRate, countTime,occThresh)
num = length(firingMap);
% num = size(firingMap,3);
infoPerSecond = zeros(num,1);
infoPerSpike = zeros(num,1);
for i  = 1:num
    if isempty(firingMap{i})
        continue;
    end
    [infoPerSecond(i), infoPerSpike(i)] = Doug_spatialInfo(firingMap{i},MeanFiringRate(i,1), countTime,occThresh);
end
% save spatialFiringInfo.mat infoPerSecond infoPerSpike;
figure
subplot(1,2,1)
bar(infoPerSecond)
xlabel('Cell #','FontSize',12,'FontName','Arial')
ylabel('infoPerSecond','FontSize',12,'FontName','Arial')
subplot(1,2,2)
bar(infoPerSpike)
xlabel('Cell #','FontSize',12,'FontName','Arial')
ylabel('infoPerSpike','FontSize',12,'FontName','Arial')
