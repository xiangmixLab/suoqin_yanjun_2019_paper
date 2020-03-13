function [place_cellsAll,place_cellsThreshAll,infoScorebootAll] = permutingSpike_parallel(file_neuronIndividuals,behavIndividuals,thresh,temp,deltaTall,occThresh,nboot,binsize,sessions,sessionUsed,neuron_lowFR,experiment)
if ~exist('occThresh','var') || isempty(occThresh)
    occThresh = 0.2;
end
if ~exist('nboot','var') || isempty(nboot)
    nboot = 100;
end
if ~exist('deltaTall','var') || isempty(deltaTall)
    deltaTall = randi([10,max(neuron.time)/1000-10],nboot,1)*1000;
end
if ~exist('temp','var') || isempty(temp)
    temp = 'S';
end
if ~exist('binsize','var') || isempty(binsize)
    binsize = 15;
end
if ~exist('sessions','var') || isempty(sessions)
    sessions = strcat('session',cellstr(num2str([1:length(behavIndividuals)])));
end
if ~exist('sessionUsed','var') || isempty(sessionUsed)
    sessionUsed = 1:length(behavIndividuals);
end
if ~exist('experiment','var') || isempty(experiment)
    experiment = 'linearTrack';
end

place_cellsAll = cell(1,length(behavIndividuals));infoScorebootAll = cell(2,length(behavIndividuals));
place_cellsThreshAll = [];
%for sessionIndex = 1:length(behavIndividuals)
for sessionIndex = 1:length(sessionUsed)
%     load([file_neuronIndividuals])
    load(file_neuronIndividuals,'neuronIndividuals')
    %neuronIndividuals = neuronIndividualsf;
    
    neuron0 = neuronIndividuals{sessionUsed(sessionIndex)};
    behav0 = behavIndividuals{sessionUsed(sessionIndex)};
    session = sessions{sessionUsed(sessionIndex)};
    [place_cells,infoScoreThresh,infoScoreSecondboot, infoScoreSpikeboot] = permutingSpike(file_neuronIndividuals,sessionUsed(sessionIndex),neuron0,behav0,thresh,temp,deltaTall,occThresh,nboot,binsize,session,neuron_lowFR,experiment);
    place_cellsAll{sessionUsed(sessionIndex)} = place_cells;
    place_cellsThreshAll(:,sessionUsed(sessionIndex)) = infoScoreThresh;
    infoScorebootAll{1,sessionUsed(sessionIndex)} = infoScoreSecondboot;
    infoScorebootAll{2,sessionUsed(sessionIndex)} = infoScoreSpikeboot;
end

