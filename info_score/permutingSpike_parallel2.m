function [place_cellsAll,place_cellsThreshAll,infoScorebootAll] = permutingSpike_parallel2(file_neuronIndividuals,behavIndividuals,thresh,temp,deltaTall,occThresh,nboot,binsize,sessions,sessionUsed,neuron_lowFR)
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

place_cellsAll = cell(1,length(behavIndividuals));infoScorebootAll = cell(2,length(behavIndividuals));
place_cellsThreshAll = [];
%for sessionIndex = 1:length(behavIndividuals)
for sessionIndex = sessionUsed
%     load([file_neuronIndividuals])
    load(file_neuronIndividuals,'neuronIndividuals')
    neuron0 = neuronIndividuals{sessionIndex};
    behav0 = behavIndividuals{sessionIndex};
    if sessionIndex >= 4
        thresh = determiningFiringEventThresh(neuron0,temp);
    end
    session = sessions{sessionIndex};
    [place_cells,infoScoreThresh,infoScoreSecondboot, infoScoreSpikeboot] = permutingSpike(file_neuronIndividuals,sessionIndex,neuron0,behav0,thresh,temp,deltaTall,occThresh,nboot,binsize,session,neuron_lowFR);
    place_cellsAll{sessionIndex} = place_cells;
    place_cellsThreshAll(:,sessionIndex) = infoScoreThresh;
    infoScorebootAll{1,sessionIndex} = infoScoreSecondboot;
    infoScorebootAll{2,sessionIndex} = infoScoreSpikeboot;
end

