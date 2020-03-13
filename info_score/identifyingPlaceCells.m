%% identify place cells
% Note: In the codes identifyingPlaceCells.m, it will load the neuronIndividuals.mat file. Please make the file name be consistent. 
addpath(genpath('E:\\Google Drive synchronization\\projects\\Project2\\Miniscope_Matlab_code\\Suoqin_CNMF-E'))
load('neuronIndividuals.mat')
load('BehavCA1_1019linearNoCNO.mat')
sectionIndex = 2;
neuron0 = neuronIndividuals{sectionIndex};
load('1019_1021combinedNeuron.mat')
thresh = determiningFiringEventThresh(neuron); %determine the neuron firing threshold
occThresh = 1; nboot = 100;
% randomly generate the dealt t for permuting the spikes
deltaTall = randi([10,890],nboot,1)*1000; % unit: ms
% deltaTall = randi([20,880],nboot,1)*1000;
[place_cells,TinfoPerSecond] = permutingSpike(sectionIndex,neuron0,behav0,thresh,deltaTall,occThresh,nboot);
save place_cells_info_10_890_CNO.mat place_cells TinfoPerSecond
