function infoScorebootAll = subsetingSpike_parallel_nooverlap(file_neuronIndividuals,behavIndividuals,thresh,temp,occThresh,nboot,binsize,sessionUsed,nbin,experiment,method_sampling)
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
    binsize=10;
end
if ~exist('sessionUsed','var') || isempty(sessionUsed)
    sessionUsed = 1:length(behavIndividuals);
end
if ~exist('nbin','var') || isempty(nbin)
    nbin = 5;
end
if ~exist('experiment','var') || isempty(experiment)
    experiment = 'linearTrack';
end

if ~exist('method_sampling','var') || isempty(method_sampling)
    method_sampling = 'jackknife';
end

infoScorebootAll = cell(2,length(behavIndividuals));
for sessionIndex = 1:length(sessionUsed)
    load(file_neuronIndividuals,'neuronIndividualsf')
    neuronIndividuals = neuronIndividualsf;
    neuron0 = neuronIndividuals{sessionUsed(sessionIndex)};
    behav0 = behavIndividuals{sessionUsed(sessionIndex)};
    %    [infoScoreSecondboot,infoScoreSpikeboot] = subsetingSpike_nooverlap_trail(file_neuronIndividuals,sessionUsed(sessionIndex),neuron0,behav0,thresh,temp,occThresh,nboot,binsize);
    if strcmpi(method_sampling,'jackknife')
        [infoScoreSecondboot,infoScoreSpikeboot] = subsetingSpike_nooverlap(file_neuronIndividuals,sessionUsed(sessionIndex),neuron0,behav0,thresh,temp,occThresh,nboot,binsize,nbin,experiment);
    elseif strcmpi(method_sampling,'bootstrap')
        [infoScoreSecondboot,infoScoreSpikeboot] = subsetingSpike_nooverlap_bootstrap(file_neuronIndividuals,sessionUsed(sessionIndex),neuron0,behav0,thresh,temp,occThresh,nboot,binsize,nbin,experiment);
    end
    infoScorebootAll{1,sessionUsed(sessionIndex)} = infoScoreSecondboot;
    infoScorebootAll{2,sessionUsed(sessionIndex)} = infoScoreSpikeboot;
end


