function [place_cells,infoScoreThresh,infoScoreSecondboot, infoScoreSpikeboot] = permutingSpike(file_neuronIndividuals,sectionIndex,neuron,behav,thresh,temp,deltaTall,occThresh,nboot,binsize,session,neuron_lowFR,experiment)
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
if ~exist('experiment','var') || isempty(experiment)
    experiment = 'linearTrack';
end

countTimeThresh = occThresh;

if strcmpi(experiment,'linearTrack')
    [firingrateAll,countAll,countTime] = calculatingCellSpatialLinearTrackForSingleData2(neuron,behav,1:size(neuron.trace,1),thresh,temp);
elseif strcmpi(experiment,'openField') | strcmpi(experiment,'OLM')
    countTimeThresh = occThresh;
    [firingrateAll,countAll,~,countTime] = calculatingCellSpatialForSingleData(neuron,behav,binsize,1:size(neuron.C,1),thresh,temp,[],[],countTimeThresh);
end

infoPerSecondnull = zeros(length(firingrateAll),1);infoPerSpikenull = infoPerSecondnull;
for j = 1:length(firingrateAll)
    MeanFiringRateAll= sum(sum(countAll{j}))/sum(sum(countTime));
    if isempty(firingrateAll{j})
        continue;
    end
    [infoPerSecondnull(j), infoPerSpikenull(j)] = Doug_spatialInfo(firingrateAll{j},MeanFiringRateAll, countTime,occThresh);
    %  [infoPerSecondnull(j), infoPerSpikenull(j)] = Doug_spatialInfo(firingrateAll{j},MeanFiringRateAll, countAll{j},occThresh);
end

% load(file_neuronIndividuals,'neuronIndividualsf')
% neuronIndividuals = neuronIndividualsf;
% load(file_neuronIndividuals,'neuronIndividuals')
%
% neuron = neuronIndividuals{sectionIndex};
% downsampling = length(neuron.time)/size(neuron.trace,2);
% neuron.time = neuron.time(1:downsampling:end);
time0 = neuron.time;
S0 = neuron.S;
trace0 = neuron.trace;
C0 = neuron.C;

infoPerSecondboot = zeros(length(firingrateAll),nboot);infoPerSpikeboot = infoPerSecondboot;
for nE = 1:nboot
    deltaT = deltaTall(nE);
    neuron.time = time0; neuron.S = S0; neuron.trace = trace0;neuron.C = C0;
    neuronboot = neuron;
    timeboot = neuronboot.time;
    %     for j = 1:size(neuron.S,1)
    %         threshI = thresh(j);
    %          idx = find(neuron.S(j,:)>threshI);
    % % idx = find(neuron.S(j,:)>0);
    %         % neuronboot.time(idx) = neuron.time(idx)+deltaT;
    %         time(idx) = time0(idx)+deltaT;
    %         index = time > max(time0);
    %         time(index) = time(index)-max(time0);
    %         [time,index] = sort(time);
    %         neuronboot.time = time;
    %         neuronboot.S(j,:) = neuronboot.S(j,index);
    %         neuronboot.trace(j,:) = neuronboot.trace(j,index);
    %     end
    
    idx = find(sum(neuron.S > repmat(thresh,1,size(neuron.S,2))));
    timeboot(idx) = time0(idx)+deltaT;
    index = timeboot > max(time0);
    timeboot(index) = timeboot(index)-max(time0);
    [timeboot,index] = sort(timeboot);
    neuronboot.time = timeboot;
    neuronboot.S = neuronboot.S(:,index);
    neuronboot.trace = neuronboot.trace(:,index);
    neuronboot.C = neuronboot.C(:,index);
    
    if strcmpi(experiment,'linearTrack')
        [firingrateAll,countAll,countTime] = calculatingCellSpatialLinearTrackForSingleData2(neuronboot,behav,1:size(neuronboot.trace,1),thresh,temp);
    elseif strcmpi(experiment,'openField') | strcmpi(experiment,'OLM')
        [firingrateAll,countAll,~,countTime] = calculatingCellSpatialForSingleData(neuronboot,behav,binsize,1:size(neuronboot.C,1),thresh,temp,[],[],countTimeThresh);
    end
    infoPerSecondbootT = zeros(length(firingrateAll),1);infoPerSpikebootT = zeros(length(firingrateAll),1);
    for j = 1:length(firingrateAll)
        MeanFiringRateAll= sum(sum(countAll{1,j}))/sum(sum(countTime));
        %         if isempty(firingrateAll{j})
        %             continue;
        %         end
        % [infoPerSecondboot(j,nE), infoPerSpikeboot(j,nE)] = Doug_spatialInfo(firingrateAll{j},MeanFiringRateAll, countTime,occThresh);
        %  [infoPerSecondboot(j,nE), infoPerSpikeboot(j,nE)] = Doug_spatialInfo(firingrateAll{j},MeanFiringRateAll, countAll{1,j},occThresh);
        [infoPerSecondbootT(j), infoPerSpikebootT(j)] = Doug_spatialInfo(firingrateAll{j},MeanFiringRateAll, countTime,occThresh);
    end
    infoPerSecondboot(:,nE) = infoPerSecondbootT; infoPerSpikeboot(:,nE) = infoPerSpikebootT;
end
infoScoreSecondboot = [infoPerSecondnull,infoPerSecondboot];
infoScoreSpikeboot = [infoPerSpikenull,infoPerSpikeboot];
% infoScore = [infoPerSecondnull,infoPerSecondboot];
% infoScore = infoPerSpikeboot;
% infoScorenull = infoPerSpikenull; 
infoScore = infoPerSecondboot;
infoScorenull = infoPerSecondnull; 

infoScore(neuron_lowFR,:) = [];

% infoScoreThresh = quantile(infoScore(:),0.95);
infoScoreThresh = quantile(infoScore,0.95,2);
% TinfoPerSecond = table([1:length(infoPerSecondnull)]',infoPerSecondnull,'VariableNames',{'neuron','infoPerSecond'});
% TinfoPerSecond2 = sortrows(TinfoPerSecond,{'infoPerSecond'},{'descend'});
% place_cells = TinfoPerSecond2.neuron(TinfoPerSecond2.infoPerSecond > infoScoreThresh);


% place_cells = TinfoPerSecond2.neuron(TinfoPerSecond2.infoPerSecond > infoScoreThresh);

TinfoScorenull = table([1:length(infoScorenull)]',infoScorenull,infoScoreThresh,'VariableNames',{'neuron','infoPerSecond','thresh'});
TinfoScorenull(neuron_lowFR,:) = [];
TinfoScorenull2 = sortrows(TinfoScorenull,{'infoPerSecond'},{'descend'});
place_cells = TinfoScorenull2.neuron(TinfoScorenull2.infoPerSecond > TinfoScorenull2.thresh);

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