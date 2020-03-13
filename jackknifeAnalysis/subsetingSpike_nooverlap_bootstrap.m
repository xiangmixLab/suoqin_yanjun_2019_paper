function [infoScoreSecondboot,infoScoreSpikeboot] = subsetingSpike_nooverlap_bootstrap(file_neuronIndividuals,sectionIndex,neuron,behav,thresh,temp,occThresh,nboot,binsize,nbin,experiment)
if ~exist('sectionIndex','var') || isempty(sectionIndex)
    sectionIndex = 1;
end
if ~exist('occThresh','var') || isempty(occThresh)
    occThresh = 0.2;
end
if ~exist('nboot','var') || isempty(nboot)
    nboot = 5;
end

if ~exist('temp','var') || isempty(temp)
    temp = 'S';
end
if ~exist('binsize','var') || isempty(binsize)
    binsize=15;
end
if ~exist('experiment','var') || isempty(experiment)
    experiment = 'linearTrack';
end
if strcmpi(experiment,'linearTrack')
    [firingrateAll,countAll,countTime] = calculatingCellSpatialLinearTrackForSingleData2(neuron,behav,1:size(neuron.trace,1),thresh,temp);
elseif strcmpi(experiment,'openField')
    countTimeThresh = occThresh;
    [firingrateAll,countAll,~,countTime] = calculatingCellSpatialForSingleData(neuron,behav,binsize,1:size(neuron.C,1),thresh,temp,[],[],countTimeThresh);
end
infoPerSecondnull = zeros(length(firingrateAll),1);infoPerSpikenull = infoPerSecondnull;
for j = 1:length(firingrateAll)
    MeanFiringRateAll= sum(sum(countAll{1,j}))/sum(sum(countTime));
    %     if isempty(firingrateAll{j})
    %         continue;
    %     end
    [infoPerSecondnull(j), infoPerSpikenull(j)] = Doug_spatialInfo(firingrateAll{j},MeanFiringRateAll, countTime,occThresh);
end

if strcmpi(experiment,'linearTrack')
    edge11 = max(behav.position(:,1))*0.05; % cut off 10% on each side
    edge12 = max(behav.position(:,1)) - edge11;
    idx1 = find(behav.position(:,1) < edge11);
    idx2 = find(behav.position(:,1) > edge12);
    idx = [idx1;idx2];
    behav.position(idx,:) = [];behav.time(idx) = [];
end

%nbin = 10;
Y = discretize(1:length(behav.time),nbin);

% load(file_neuronIndividuals,'neuronIndividualsf')
% neuronIndividuals = neuronIndividualsf;
% neuron = neuronIndividuals{sectionIndex};
time0 = neuron.time;
S0 = neuron.S;
trace0 = neuron.trace;
C0 = neuron.C;

neuron.pos = interp1(behav.time,behav.position,neuron.time);
pos0 = neuron.pos;

infoPerSecondboot = zeros(length(firingrateAll),nboot);infoPerSpikeboot = infoPerSecondboot;
for nE = 1:nboot
    %     downsampling = length(neuron.time)/size(neuron.trace,2);
    %     neuron.time = neuron.time(1:downsampling:end);
    
    neuron.time = time0; neuron.S = S0; neuron.trace = trace0; neuron.C = C0;neuron.pos = pos0;
    
    neuronboot = neuron;
    behavboot = behav;
    
    bins_use = sort(randi(nbin,1,nbin));
    bins_use_sample = bins_use;
    frame_pass = [];
    for i = 1:length(bins_use_sample)
        frame_pass = [frame_pass,find(Y == bins_use_sample(i))];
    end
    
    % find the nearest points based on the time
    %     idx_neuron = knnsearch(neuron.time,behav.time([frame_pass(1),frame_pass(end)]));
    %     frame_pass_neuron = idx_neuron(1):idx_neuron(end);
    idx_neuron = knnsearch(neuron.time,behav.time);
    frame_pass_neuron = idx_neuron;
    
    behavboot.position = behavboot.position(frame_pass,:);
    behavboot.time = behavboot.time(frame_pass);
    
    neuronboot.S = neuronboot.S(:,frame_pass_neuron);
    neuronboot.trace = neuronboot.trace(:,frame_pass_neuron);
    neuronboot.C = neuronboot.C(:,frame_pass_neuron);
    neuronboot.time = neuronboot.time(frame_pass_neuron);
    neuronboot.num2read = length(frame_pass_neuron);
    neuronboot.pos = neuronboot.pos(frame_pass_neuron,:);
    if strcmpi(experiment,'linearTrack')
        [firingrateAll,countAll,countTime] = calculatingCellSpatialLinearTrackForSingleData3(neuronboot,behavboot,1:size(neuronboot.trace,1),thresh,temp);
    elseif strcmpi(experiment,'openField')
        [firingrateAll,countAll,~,countTime] = calculatingCellSpatialForSingleData(neuronboot,behavboot,binsize,1:size(neuronboot.C,1),thresh,temp,[],[],countTimeThresh);
    end
    infoPerSecondbootT = zeros(length(firingrateAll),1);infoPerSpikebootT = zeros(length(firingrateAll),1);
    for j = 1:length(firingrateAll)
        MeanFiringRateAll= sum(sum(countAll{1,j}))/sum(sum(countTime));
        %         if isempty(firingrateAll{j})
        %             continue;
        %         end
        %         [infoPerSecondboot(j,nE), infoPerSpikeboot(j,nE)] = Doug_spatialInfo(firingrateAll{j},MeanFiringRateAll, countTime,occThresh);
        [infoPerSecondbootT(j), infoPerSpikebootT(j)] = Doug_spatialInfo(firingrateAll{j},MeanFiringRateAll, countTime,occThresh);
    end
    infoPerSecondboot(:,nE) = infoPerSecondbootT; infoPerSpikeboot(:,nE) = infoPerSpikebootT;
end

infoScoreSecondboot = [infoPerSecondnull,infoPerSecondboot];

infoScoreSpikeboot = [infoPerSpikenull,infoPerSpikeboot];
