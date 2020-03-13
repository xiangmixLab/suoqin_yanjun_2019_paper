    %%% temporal dynamical analysis of Ca2+ imaging data
addpath(genpath('dynamicalAnalysiCodes'))
%% load the combined neuron data file
load('1019_1021combinedNeuron.mat')
%% given a threshold for each neuron.firing
thresh = determiningFiringEventThresh(neuron,'trace'); % 10% of the peak
%% load the individual neuron data and extract the time information in each condition
load('neuronIndividuals.mat')
neuron1 = neuronIndividuals{1};
neuron2 = neuronIndividuals{2};
behav1 = importdata('BehavCA1_1019linearNoCNO.mat');
behav2 = importdata('BehavCA1_1021_linearCNO.mat');
% [activityEnsem1,activityEnsemLR1,activityEnsemRL1] = calculatingEnsembleActivityLinearTrack(neuron1,behav1,'trace');
% [activityEnsem2,activityEnsemLR2,activityEnsemRL2] = calculatingEnsembleActivityLinearTrack(neuron2,behav2,'trace');
% %

[ensemIfo1,ensemIfoLR1,ensemIfoRL1] = calculatingEnsembleActivityLinearTrack(neuron1,behav1,thresh,'S');
[ensemIfo2,ensemIfoLR2,ensemIfoRL2] = calculatingEnsembleActivityLinearTrack(neuron2,behav2,thresh,'S');

ensemIfoLR1.FR(isnan(ensemIfoLR1.FR)) = 0;ensemIfoRL1.FR(isnan(ensemIfoRL1.FR)) = 0;
d = pdist2(ensemIfoLR2.FR',ensemIfoRL2.FR','correlation');
% d = pdist2(ensemIfo1.FR',ensemIfo2.FR','correlation');
corr = 1-d;
% corr = diag(corr);
figure;
imagesc((corr))
colormap(flipud(hot))
c = colorbar;

d = pdist2(ensemIfo1.count',ensemIfo2.count','correlation');
corr = 1-d;
figure;
imagesc(corr)
colormap(flipud(hot))
c = colorbar;

d = pdist2(ensemIfoLR1.activity',ensemIfoLR2.activity','correlation');
corr = 1-d;
figure;
imagesc(corr)
colormap(flipud(hot))
c = colorbar;

%% 
[ensemIfo1,ensemIfoLR1,ensemIfoRL1,ensemIfoLRodd1,ensemIfoRLodd1,ensemIfoLReven1,ensemIfoRLeven1] = ...
    calculatingEnsembleActivityLinearTrack2(neuron1,behav1,thresh,'trace');
[ensemIfo2,ensemIfoLR2,ensemIfoRL2,ensemIfoLRodd2,ensemIfoRLodd2,ensemIfoLReven2,ensemIfoRLeven2] = ...
    calculatingEnsembleActivityLinearTrack2(neuron2,behav2,thresh,'trace');

d = pdist2(ensemIfoRLodd1.FR',ensemIfoRLeven1.FR','correlation');
corr = 1-d;
mean(diag(corr))
figure;
imagesc((corr))
colormap(flipud(hot))
c = colorbar;
xlabel('Even trails','FontSize',10);ylabel('Odd trails','FontSize',10)
title('Correlation map based on FR (CNO)','FontSize',10)
axis square

d = pdist2(ensemIfoRLodd2.activity',ensemIfoRLeven2.activity','correlation');
corr = 1-d;
mean(diag(corr))
figure;
imagesc((corr))
colormap(flipud(hot))
c = colorbar;
xlabel('Even trails','FontSize',10);ylabel('Odd trails','FontSize',10)
title('Correlation map based on Activity (CNO)','FontSize',10)
axis square




