%% after the neuron and behav data have been processed, continue with this process to get the spatial coding map of the neurons
% This process is for individual sessions, not combined sessions with
% different experiment conditions
addpath(genpath('E:\Google Drive synchronization\projects\Project2\Miniscope_Matlab_code'))
%% 1st step: load the varibale neuron and behav
%% 2nd step
thresh = determiningFiringEventThresh(neuron,'S'); %determine the neuron firing threshold
%% plot all the trace (optional)
plottingTrace(neuron,1:size(neuron.trace,1),1)
%% 3rd step, incorporate the firing rate information and then plot the final spatial coding heat map
mice = 'M3244F';sessions = {'0217Ctrl','0221CNO','0224PCtrl'};
sessionIndex = 1;
neuron0 = neuronIndividuals{sessionIndex};
behav0 = behavIndividuals{sessionIndex};
% [firingrateAll,countAll,countTime] = calculatingCellSpatialForSingleData(neuron,behav,1:size(neuron.trace,1),thresh);
[firingrateAll,countAll,countTime] = calculatingCellSpatialLinearTrackForSingleData2(neuron0,behav0,1:size(neuron0.trace,1),thresh,'S');

% plottingFiringBehaviorSpatialForSingleData(neuron0,behav0,firingrateAll,1:size(neuron.trace,1),thresh,'S',20,[mice,sessions{sessionIndex}])
plottingFiringBehaviorSpatialForSingleData(neuron0,behav0,firingrateAll,1:size(neuron.trace,1),thresh,'S',5,[mice,sessions{sessionIndex}])
%% Calculate mean firing rate for all neurons
MeanFiringRateAll=zeros(size(firingrateAll,2),1);
for i=1:size(firingrateAll,2)
    MeanFiringRateAll(i,1)= sum(sum(countAll{1,i}))/sum(sum(countTime));
end
%% Calculate the info score
[infoPerSecond, infoPerSpike] = comparisonSpatialInfo(firingrateAll, MeanFiringRateAll, countTime,1);
save([mice,'_',sessions{sessionIndex},'_','spatialFiringInfo.mat'],'firingrateAll','countAll','countTime')
save([mice,'_',sessions{sessionIndex},'_','MeanFiringRateAll.mat'],'MeanFiringRateAll')
save([mice,'_',sessions{sessionIndex},'_','InfoPerSecAll.mat'],'infoPerSecond')
save([mice,'_',sessions{sessionIndex},'_','InfoPerSpikeAll.mat'],'infoPerSpike')

%% compare the info scores across sessions
info_scoreAll = [1:length(infoPerSecond)]';
for i = 1:length(sessions)
    load([mice,'_',sessions{i},'_','InfoPerSecAll.mat'])
    info_scoreAll = [info_scoreAll,infoPerSecond];
end
[~,order] = sort(info_scoreAll(:,2),'descend');
info_scoreAll = info_scoreAll(order,:);
save info_score_comparisons_0217Ctrl_0221CNO_0224PCtrl.mat info_scoreAll
figure
imagesc(info_scoreAll(:,2:end))
colormap(flipud(hot))
xlabels = {'Contrl2','CNO1', 'Post-contrl1'};
set(gca,'Xtick',1:length(xlabels))
set(gca,'Xticklabel',xlabels,'FontSize',8);
c = colorbar;
c.Location = 'east';
c.Label.String = 'Info Per Second';
c.Label.FontSize = 8;%c.Label.FontWeight = 'bold';
c.FontSize = 7;
c.Position = [0.7 .11 .02 .2];
%% identify place cells
% Note: In the codes identifyingPlaceCells.m, it will load the neuronIndividuals.mat file. Please make the file name be consistent.
% load('neuronIndividuals.mat')
% load('BehavCA1_1019linearNoCNO.mat')
% sessionIndex = 1;
% neuron0 = neuronIndividuals{sessionIndex};
% thresh = determiningFiringEventThresh(neuron0); %determine the neuron firing threshold
occThresh = 1; nboot = 100;
% randomly generate the dealt t for perpute the spike
deltaTall = randi([10,890],nboot,1)*1000; % unit: ms
% deltaTall = randi([20,880],nboot,1)*1000;
%% Note1: one should change the code in line #29 in permutingSpike, in order to load the right neuronIndividual file
% load('neuronIndividuals_filtered_M3244F_0217Ctrl_0221CNO_0224PCtrl.mat')
%% Note 2: After running the permutingSpike, please reload the neurionIndividual variable if you are going to analyze the next session.
%%
[place_cells,TinfoPerSecond] = permutingSpike(sessionIndex,neuron0,behav0,thresh,'S',deltaTall,occThresh,nboot);
save([mice,'_',sessions{sessionIndex},'_','place_cells_info_10_890.mat'],'place_cells', 'TinfoPerSecond')

%% correlation of individual neuron
% load data
mice = 'M3244F';sessions = {'0217Ctrl','0221CNO','0224PCtrl'};
FR_Ifo = cell(1,3);
for i = 1:3
    load([mice,'_',sessions{i},'_','spatialFiringInfo.mat'],'firingrateAll');
    FR_Ifo{i} = firingrateAll;
end

FR_corr = zeros(length(FR_Ifo{1}),3);
for k = 1:length(FR_Ifo{1})
    kk = 0;
    for i = 1:length(FR_Ifo)-1
        FR_Ifo{i}{k}(isnan(FR_Ifo{i}{k})) = 0;
        for j = i+1:length(FR_Ifo)
            kk = kk + 1;
            FR_Ifo{j}{k}(isnan(FR_Ifo{j}{k})) = 0;
            FR_corr(k,kk) = corr2(FR_Ifo{i}{k}(3:4,end-27:end),FR_Ifo{j}{k}(3:4,end-27:end));
        end
    end
end
load('M3244F_0217Ctrl_place_cells_info_10_890.mat')
figure
imagesc(FR_corr(TinfoPerSecond.neuron,:))
colormap(flipud(hot))
xlabels = {'Ctrl-CNO','Ctrl-PCtrl', 'CNO-PCtrl'};
set(gca,'Xtick',1:length(xlabels))
set(gca,'Xticklabel',xlabels,'FontSize',8);
% xtickangle(45)
% set(gca,'Ytick',1:length(xlabels))
% set(gca,'Yticklabel',xlabels,'FontSize',8);
% ytickangle(45)
% c = colorbar;
% c.Location = 'east';
% c.Label.String = 'Correlation';
% c.Label.FontSize = 8;%c.Label.FontWeight = 'bold';
% c.FontSize = 7;
% c.Position = [0.7 .11 .02 .2];

% distribution of correlations
FR_corr_ordered = FR_corr(TinfoPerSecond.neuron,:);
figure
clf
subplot(1,2,1)
% clf
% boxplot(FR_corr_ordered(1:30,:))
h = violinplot(FR_corr_ordered(1:50,:));
xlim([0.5 length(h)+0.5])
%     for j = 1:length(h)
%         h(j).ViolinColor = colorTSNE(j,:);
%     end
xlabels = {'Ctrl-CNO','Ctrl-PCtrl', 'CNO-PCtrl'};
set(gca,'Xtick',1:length(xlabels))
set(gca,'Xticklabel',xlabels,'FontSize',8);
ylabel('Spatial correlation','FontSize',10)
subplot(1,2,2)
% boxplot(FR_corr_ordered(end-29:end,:))
h = violinplot(FR_corr_ordered(end-49:end,:));
xlim([0.5 length(h)+0.5])
set(gca,'Xtick',1:length(xlabels))
set(gca,'Xticklabel',xlabels,'FontSize',8);
ylabel('Spatial correlation','FontSize',10)

