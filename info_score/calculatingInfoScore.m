%% after the neuron and behav data have been processed, continue with this process to get the spatial coding map of the neurons
% This process is for individual sessions, not combined sessions with different experiment conditions
%% 1st step: load the varibale neuron and behav
%% 2nd step
thresh = determiningFiringEventThresh(neuron); %determine the neuron firing threshold
%% plot all the trace (optional)
% plottingTrace(neuron,1:size(neuron.trace,1),0)

%% 3rd step, incorporate the firing rate information and then plot the final spatial coding heat map
[firingrateAll,countAll,countTime] = calculatingCellSpatialForSingleData(neuron,behav,1:size(neuron.trace,1),thresh);
% [firingrateAll,countAll,countTime] = calculatingCellSpatialLinearTrackForSingleData(neuron,behav,1:size(neuron.trace,1),thresh);

plottingFiringBehaviorSpatialForSingleData(neuron,behav,firingrateAll,1:size(neuron.trace,1),thresh,5,'M3321_Ctrl1')

%% Calculate mean firing rate for all neurons
MeanFiringRateAll=zeros(size(firingrateAll,2),1);
for i=1:size(firingrateAll,2)
    MeanFiringRateAll(i,1)= sum(sum(countAll{1,i}))/sum(sum(countTime));
end

%% Calculate the info score
[infoPerSecond, infoPerSpike] = comparisonSpatialInfo(firingrateAll, MeanFiringRateAll, countTime,1);

save spatialInfoScore.mat infoPerSecond infoPerSpike;
