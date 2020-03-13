function [FiringRateAll,CountAll,CountTimeAll,CountTimeAll1,MeanFiringRateAll,InfoScoreAll,InfoSpikeAll,AmplitudeAll,binInfoAll] = calculating_FR_IS_parallel(neuronIndividuals,behavIndividuals,binsize,temp,thresh,countTimeThresh,experiment)
if ~exist('temp','var') || isempty(temp)
    temp = 'S';
end
if ~exist('experiment','var') || isempty(experiment)
    experiment = 'linearTrack';
end
num = size(neuronIndividuals{1}.C,1);
InfoScoreAll = [1:num]';InfoSpikeAll = [1:num]';
MeanFiringRateAll = [1:num]';
FiringRateAll = cell(num,length(neuronIndividuals));CountAll = cell(num,length(neuronIndividuals));
CountTimeAll = cell(num,length(neuronIndividuals));
CountTimeAll1 = cell(1,length(neuronIndividuals));
binInfoAll = cell(1,length(neuronIndividuals));
AmplitudeAll = cell(num,length(neuronIndividuals));
for sessionIndex = 1:length(neuronIndividuals)
    neuron0 = neuronIndividuals{sessionIndex};
    behav0 = behavIndividuals{sessionIndex};
    
    if strcmpi(experiment,'linearTrack')
        [firingRate,count,countTime] = calculatingCellSpatialLinearTrackForSingleData2(neuron0,behav0,1:size(neuron0.C,1),thresh,temp);
    elseif strcmpi(experiment,'openField')
        countTimeThresh = occThresh;
        [firingRate,count,~, countTime] = calculatingCellSpatialForSingleData(neuron0,behav0,binsize,1:size(neuron0.C,1),thresh,temp,[],[],countTimeThresh); %OPEN FIELD
    elseif strcmpi(experiment,'OLM')
        [firingRate,count,countTime2, countTime,amplitude,binInfo] = calculatingCellSpatialForSingleData(neuron0,behav0,binsize,1:size(neuron0.C,1),thresh,temp,[],[],countTimeThresh); % OLM
    end
    
end

FiringRateAll(:,sessionIndex) = firingRate';
CountAll(:,sessionIndex) = count';
CountTimeAll1{sessionIndex} = countTime;

if strcmpi(experiment,'OLM')
    CountTimeAll(:,sessionIndex) = countTime2';
    AmplitudeAll(:,sessionIndex) = amplitude';
    binInfoAll{sessionIndex} = binInfo;
end
%% Calculate mean firing rate for all neurons
MeanFiringRate=zeros(size(firingRate,2),1);
for i=1:size(firingRate,2)
    MeanFiringRate(i,1)= sum(sum(count{1,i}))/sum(sum(countTime));
end
%% Calculate the info score
[infoPerSecond, infoPerSpike] = comparisonSpatialInfo(firingRate, MeanFiringRate, countTime,countTimeThresh(1));
%% combine the info scores across sessions
InfoScoreAll = [InfoScoreAll,infoPerSecond];
InfoSpikeAll = [InfoSpikeAll,infoPerSpike];
MeanFiringRateAll = [MeanFiringRateAll,MeanFiringRate];

% [~,order] = sort(info_scoreAll(:,2),'descend');
% info_scoreAll = info_scoreAll(order,:);
% MeanFiringRateAll = MeanFiringRateAll(order,:);

