function [infoScoreSecondboot,infoScoreSpikeboot] = subsetingSpike_nooverlap_trail(file_neuronIndividuals,sectionIndex,neuron,behav,thresh,temp,occThresh,nboot,binsize)
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

[firingrateAll,countAll,countTime] = calculatingCellSpatialLinearTrackForSingleData2(neuron,behav,1:size(neuron.trace,1),thresh,temp);
%[firingrateAll,countAll,countTimeAll,countTime] = calculatingCellSpatialForSingleData(neuron,behav,binsize,1:size(neuron.C,1),thresh,temp,[],[],countTimeThresh);
infoPerSecondnull = zeros(length(firingrateAll),1);infoPerSpikenull = infoPerSecondnull;
for j = 1:length(firingrateAll)
    MeanFiringRateAll= sum(sum(countAll{1,j}))/sum(sum(countTime));
    if isempty(firingrateAll{j})
        continue;
    end
    [infoPerSecondnull(j), infoPerSpikenull(j)] = Doug_spatialInfo(firingrateAll{j},MeanFiringRateAll, countTime,occThresh);
end

[Y,E] = discretize(1:max(neuron.trialNum),nboot);E(1) = 1;
trialNum_bins = zeros(length(neuron.trialNum),1);
for i = 1:nboot
    idx1 = neuron.trialNum >= E(i);
    idx2 = neuron.trialNum < E(i+1);
    idx = idx1 & idx2;
    trialNum_bins(idx) = i;
end
% deltaTall = randi([10,890],nboot,1)*1000;
% deltaTall = randi([20,880],nboot,1)*1000;
load(file_neuronIndividuals,'neuronIndividuals')
neuron = neuronIndividuals{sectionIndex};
time0 = neuron.time;
S0 = neuron.S;
trace0 = neuron.trace;
C0 = neuron.C;
% 
neuron.pos = interp1(behav.time,behav.position,neuron.time);
pos0 = neuron.pos;
infoPerSecondboot = zeros(length(firingrateAll),nboot);infoPerSpikeboot = infoPerSecondboot;
for nE = 1:nboot
    
    %     downsampling = length(neuron.time)/size(neuron.trace,2);
    %     neuron.time = neuron.time(1:downsampling:end);
    neuron.time = time0; neuron.S = S0; neuron.trace = trace0; neuron.C = C0;neuron.pos = pos0; 
    neuronboot = neuron;
    behavboot = behav;
    
    % frame_pass_neuron = find(trialNum_bins == nE);
    bins_use = setdiff(1:nboot,nE);
%     bins_use_sample = bins_use(randperm(length(bins_use)));
   bins_use_sample = bins_use;
    frame_pass_neuron = [];
    for i = 1:length(bins_use_sample)
        frame_pass_neuron = [frame_pass_neuron;find(trialNum_bins == bins_use_sample(i))];
    end
    % find the nearest points based on the time
    %     idx_neuron = knnsearch(neuron.time,behav.time([frame_pass(1),frame_pass(end)]));
    %     frame_pass_neuron = idx_neuron(1):idx_neuron(end);
    frame_pass = knnsearch(behav.time,neuron.time(frame_pass_neuron));
    
    behavboot.position = behavboot.position(frame_pass,:);
    behavboot.time = behavboot.time(frame_pass);
    
    neuronboot.S = neuronboot.S(:,frame_pass_neuron);
    neuronboot.trace = neuronboot.trace(:,frame_pass_neuron);
    neuronboot.C = neuronboot.C(:,frame_pass_neuron);
    neuronboot.time = neuronboot.time(frame_pass_neuron);
    neuronboot.num2read = length(frame_pass_neuron);
    neuronboot.pos = neuronboot.pos(frame_pass_neuron,:);
    
    [firingrateAll,countAll,countTime] = calculatingCellSpatialLinearTrackForSingleData3(neuronboot,behavboot,1:size(neuronboot.trace,1),thresh,temp);
    %  [firingrateAll,countAll,countTimeAll,countTime] = calculatingCellSpatialForSingleData(neuronboot,behavboot,binsize,1:size(neuronboot.C,1),thresh,temp,[],[],countTimeThresh);
    for j = 1:length(firingrateAll)
        MeanFiringRateAll= sum(sum(countAll{1,j}))/sum(sum(countTime));
        if isempty(firingrateAll{j})
            continue;
        end
        [infoPerSecondboot(j,nE), infoPerSpikeboot(j,nE)] = Doug_spatialInfo(firingrateAll{j},MeanFiringRateAll, countTime,occThresh);
    end
end
infoScoreSecondboot = [infoPerSecondnull,infoPerSecondboot];

infoScoreSpikeboot = [infoPerSpikenull,infoPerSpikeboot];

% infoScore = [infoPerSecondnull,infoPerSecondboot];
% infoScoreThresh = quantile(infoScore(:),0.95);
% infoScore = [infoPerSecondnull,infoPerSecondboot];

% place_cells = find(infoPerSecondnull > infoScoreThresh);

% TinfoPerSecond = table([1:length(infoPerSecondnull)]',infoPerSecondnull,'VariableNames',{'neuron','infoPerSecond'});
% TinfoPerSecond2 = sortrows(TinfoPerSecond,{'infoPerSecond'},{'descend'});
% place_cells = TinfoPerSecond2.neuron(TinfoPerSecond2.infoPerSecond > infoScoreThresh);

% folderName = fullfile('results','figures');
% if ~exist(folderName, 'dir')
%     mkdir(folderName);
% end
%
% hFig = figure('position', [200, 200, 220,220]);
% histogram(infoPerSecondboot(:),'Normalization','probability');
% hold on
% ylim = get(gca,'ylim');
% line([infoScoreThresh infoScoreThresh],get(gca,'ylim'),'LineStyle','--','Color','r','LineWidth',1)
% xlim([min(infoScore(:)),max(infoScore(:))])
% text(infoScoreThresh+0.05,range(ylim)/2,['Thresh = ',num2str(infoScoreThresh,'%.4f')],'Color','red','FontSize',8)
% xlabel('Spatial info score','FontSize',10)
% ylabel('Frequency','FontSize',10)
% title([session,': n=',num2str(length(infoPerSecondnull)),', ', '#place cell: ', num2str(length(place_cells))],'FontSize',10)
% saveas(gcf,fullfile(folderName,['shuffled_distribution_info_score_',session,'.pdf']))
% saveas(gcf,fullfile(folderName,['shuffled_distribution_info_score_',session,'.fig']))
