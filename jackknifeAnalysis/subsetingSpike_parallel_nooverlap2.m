function infoScorebootAll = subsetingSpike_parallel_nooverlap2(file_neuronIndividuals,behavIndividuals,thresh,temp,occThresh,nboot,binsize,sessionUsed)
if ~exist('occThresh','var') || isempty(occThresh)
    occThresh = 0.2;
end
if ~exist('nboot','var') || isempty(nboot)
    nboot = 100;
end

if ~exist('temp','var') || isempty(temp)
    temp = 'S';
end
if ~exist('binsize','var') || isempty(binsize)
    binsize=15;
end

if ~exist('sessionUsed','var') || isempty(sessionUsed)
    sessionUsed = 1:length(behavIndividuals);
end

infoScorebootAll = cell(2,length(behavIndividuals));
for sessionIndex = 1:length(sessionUsed)
%for sessionIndex = 1
    %     load([file_neuronIndividuals])
    load(file_neuronIndividuals,'neuronIndividualsf')
    neuronIndividuals = neuronIndividualsf;
    neuron0 = neuronIndividuals{sessionUsed(sessionIndex)};
    behav0 = behavIndividuals{sessionUsed(sessionIndex)};
    if sessionUsed(sessionIndex) >= 4
        thresh = determiningFiringEventThresh(neuron0,temp);
    end
%    [infoScoreSecondboot,infoScoreSpikeboot] = subsetingSpike_nooverlap_trail(file_neuronIndividuals,sessionUsed(sessionIndex),neuron0,behav0,thresh,temp,occThresh,nboot,binsize);
     [infoScoreSecondboot,infoScoreSpikeboot] = subsetingSpike_nooverlap(file_neuronIndividuals,sessionUsed(sessionIndex),neuron0,behav0,thresh,temp,occThresh,nboot,binsize);     
 %[infoScoreSecondboot,infoScoreSpikeboot] = subsetingSpike_nooverlap_bootstrap(file_neuronIndividuals,sessionUsed(sessionIndex),neuron0,behav0,thresh,temp,occThresh,nboot,binsize);     
 infoScorebootAll{1,sessionIndex} = infoScoreSecondboot;
    infoScorebootAll{2,sessionIndex} = infoScoreSpikeboot;
end

