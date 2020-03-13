function [place_cells,infoScoreThresh,infoScoreSecondboot, infoScoreSpikeboot] = permutingSpike_individual(file_neuronIndividuals,sectionIndex,neuron,behav,thresh,temp,deltaTall,occThresh,nboot,binsize,session,neuron_lowFR)
if ~exist('sectionIndex','var') || isempty(sectionIndex)
    sectionIndex = 1;
end
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
    temp = 'trace';
end
if ~exist('binsize','var') || isempty(binsize)
    binsize=15;
end
if ~exist('session','var') || isempty(session)
    session = ' ';
end
if ~exist('neuron_lowFR','var') || isempty(neuron_lowFR)
    neuron_lowFR = [];
end

countTimeThresh = occThresh;

[firingrateAll,countAll,countTime] = calculatingCellSpatialLinearTrackForSingleData2(neuron,behav,1:size(neuron.trace,1),thresh,temp);
%[firingrateAll,countAll,countTimeAll,countTime] = calculatingCellSpatialForSingleData(neuron,behav,binsize,1:size(neuron.C,1),thresh,temp,[],[],countTimeThresh);
infoPerSecondnull = zeros(length(firingrateAll),1);infoPerSpikenull = infoPerSecondnull;
for j = 1:length(firingrateAll)
    MeanFiringRateAll= sum(sum(countAll{j}))/sum(sum(countTime));
    if isempty(firingrateAll{j})
        continue;
    end
    [infoPerSecondnull(j), infoPerSpikenull(j)] = Doug_spatialInfo(firingrateAll{j},MeanFiringRateAll, countTime,occThresh);
    %  [infoPerSecondnull(j), infoPerSpikenull(j)] = Doug_spatialInfo(firingrateAll{j},MeanFiringRateAll, countAll{j},occThresh);
end
num = length(firingrateAll);
infoPerSecondboot = zeros(num,nboot);infoPerSpikeboot = infoPerSecondboot;
for nE = 1:nboot
    deltaT = deltaTall(nE);
    load(file_neuronIndividuals,'neuronIndividuals')
    neuron = neuronIndividuals{sectionIndex};
    downsampling = length(neuron.time)/size(neuron.trace,2);
    neuron.time = neuron.time(1:downsampling:end);
    time0 = neuron.time;
    neuronboot = neuron;
    dataS = neuron.S;
    infoPerSecondbootT = zeros(num,1);infoPerSpikebootT = zeros(num,1);
    for j = 1:size(dataS,1)
        timeboot = time0;
        threshI = thresh(j);
        idx = find(dataS(j,:)>threshI);
        timeboot(idx) = time0(idx)+deltaT;
        index = timeboot > max(time0);
        timeboot(index) = timeboot(index)-max(time0);
        [timeboot,index] = sort(timeboot);
        neuronboot.time = timeboot;
        neuronboot.S(j,:) = neuronboot.S(j,index);
        neuronboot.trace(j,:) = neuronboot.trace(j,index);
        [firingrateAll,countAll,countTime] = calculatingCellSpatialLinearTrackForSingleData2(neuronboot,behav,j,thresh,temp);
        MeanFiringRateAll= sum(sum(countAll{1}))/sum(sum(countTime));
        [infoPerSecondbootT(j), infoPerSpikebootT(j)] = Doug_spatialInfo(firingrateAll{1},MeanFiringRateAll, countTime,occThresh);
    end
    
    infoPerSecondboot(:,nE) = infoPerSecondbootT; infoPerSpikeboot(:,nE) = infoPerSpikebootT;
end
infoScoreSecondboot = [infoPerSecondnull,infoPerSecondboot];
infoScoreSpikeboot = [infoPerSpikenull,infoPerSpikeboot];
% infoScore = [infoPerSecondnull,infoPerSecondboot];
infoScore = infoPerSpikeboot;

infoScore(neuron_lowFR,:) = [];

% infoScoreThresh = quantile(infoScore(:),0.95);
infoScoreThresh = quantile(infoScore,0.95,2);
% TinfoPerSecond = table([1:length(infoPerSecondnull)]',infoPerSecondnull,'VariableNames',{'neuron','infoPerSecond'});
% TinfoPerSecond2 = sortrows(TinfoPerSecond,{'infoPerSecond'},{'descend'});
% place_cells = TinfoPerSecond2.neuron(TinfoPerSecond2.infoPerSecond > infoScoreThresh);

TinfoPerSecond = table([1:length(infoPerSpikenull)]',infoPerSpikenull,infoScoreThresh,'VariableNames',{'neuron','infoPerSecond','thresh'});
TinfoPerSecond(neuron_lowFR,:) = [];
TinfoPerSecond2 = sortrows(TinfoPerSecond,{'infoPerSecond'},{'descend'});
% place_cells = TinfoPerSecond2.neuron(TinfoPerSecond2.infoPerSecond > infoScoreThresh);
place_cells = TinfoPerSecond2.neuron(TinfoPerSecond2.infoPerSecond > TinfoPerSecond2.thresh);

%{
folderName = fullfile('results','figures');
if ~exist(folderName, 'dir')
    mkdir(folderName);
end

hFig = figure('position', [200, 200, 220,220]);
histogram(infoScore(:),'Normalization','probability');
hold on
ylim = get(gca,'ylim');
line([infoScoreThresh infoScoreThresh],get(gca,'ylim'),'LineStyle','--','Color','r','LineWidth',1)
xlim([min(infoScore(:)),max(infoScore(:))])
text(infoScoreThresh+0.05,range(ylim)/2,['Thresh = ',num2str(infoScoreThresh,'%.4f')],'Color','red','FontSize',8)
xlabel('Spatial info score','FontSize',10)
ylabel('Frequency','FontSize',10)
title([session,': n=',num2str(length(infoPerSecondnull)),', ', '#place cell: ', num2str(length(place_cells))],'FontSize',10)
saveas(gcf,fullfile(folderName,['shuffled_distribution_info_score_',session,'.pdf']))
saveas(gcf,fullfile(folderName,['shuffled_distribution_info_score_',session,'.fig']))
%}