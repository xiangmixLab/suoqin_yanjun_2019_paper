%%%%%%%%% Calcium imaging data analysis (Jan-2019)-(updated April-2019)-(updated July-2019) %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% README start %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is a general pipeline for performing jackknife or bootstrap sampling analysis, including linear track and open field
% For open field, we used different parameters such as binsize, occThresh and different methods for determining the neuron firing threshold. 
% This code load neuron and behavior variable files from 'filefolder' and save results into 'filefolder'
% Functionality of this code: recalculate info score using jackknife or bootstrap sampling method
%%%%%%%%%%%%%%%%%%%% README end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [infoScorebootAll, sessionsSelected] = demo_calciumImagingData_analysis_openField_subseting(filefolder,mice,method_sampling)
if ~exist('method_sampling','var') || isempty(method_sampling)
    method_sampling = 'jackknife';
end

folderName = fullfile(filefolder,'results','figures');
if ~exist(folderName, 'dir')
    mkdir(folderName);
end
cd(filefolder)

%% open field
% load individual neuron file and behavior file
file_neuronIndividuals = fullfile(filefolder,['neuronIndividualsf.mat'])
load(file_neuronIndividuals)
neuronIndividuals = neuronIndividualsf;
load(fullfile(filefolder,['behavIndividualsf.mat']))
behavIndividuals = behavIndividualsf;
sessionsSelected = [1 2 3];
sessionCtrl = 1;
sessions = strcat('session',cellstr(num2str([1:3]')));

% load combined neuron variable for determining the firing threshold
neuronfile = fullfile(filefolder,[mice,'_Neuron_filter.mat'])
load(neuronfile)

%% ------------ analysis start: ------------------
%% determine the neuron firing threshold
temp = 'S';
thresh = determiningFiringEventThresh(neuron,temp,'sd');

%% calculate the firing rate and information score information
%% incorporate the firing rate information and then plot the final spatial coding heat map
binsize = 2; countTimeThresh = 0.1;% open filed
% [FiringRateAll,CountAll,~,CountTimeAll1,MeanFiringRateAll,InfoScoreAll,InfoSpikeAll] = calculating_FR_IS_parallel(neuronIndividuals,behavIndividuals,binsize,temp,thresh,countTimeThresh);
% save(fullfile(filefolder,'results',['FiringRate_',mice,'.mat']), 'FiringRateAll','CountAll','CountTimeAll1','MeanFiringRateAll')
% save(fullfile(filefolder,'results',['InfoScore_',mice,'.mat']), 'InfoScoreAll','InfoSpikeAll')

%% identify place cells
% % occThresh = 1; nboot = 100;
% occThresh = 0.2; nboot = 1000;
% % randomly generate the dealt t for perpute the spike
% deltaTall = randi([10,890],nboot,1)*1000; % unit: ms
% % deltaTall = randi([20,880],nboot,1)*1000;
% % neuron_lowFR = find(sum(MeanFiringRateAll(:,2:end) < 0.1,2) >= 2);
% neuron_lowFR = [];
%
% [place_cellsAll,place_cellsThreshAll,infoScorebootAll] = permutingSpike_parallel(file_neuronIndividuals,behavIndividuals,thresh,temp,deltaTall,occThresh,nboot,binsize,sessions,sessionCtrl,neuron_lowFR);
% % [place_cells,TinfoPerSecond] = permutingSpike(file_neuronIndividuals,sessionIndex,neuron0,behav0,thresh,'S',deltaTall,occThresh,nboot,binsize);
% save(fullfile(filefolder,'results',['place_cells_',mice,'.mat']), 'place_cellsAll','place_cellsThreshAll','infoScorebootAll')

occThresh = 0.1; nbin = 10; nboot = 1000; experiment = 'openField';
infoScorebootAll = subsetingSpike_parallel_nooverlap(file_neuronIndividuals,behavIndividuals,thresh,temp,occThresh,nboot,binsize,sessionsSelected,nbin,experiment, method_sampling);
save(fullfile(filefolder,['infoScorebootAll_subsetingBiningTime_',num2str(nbin),'_',num2str(nboot),mice,'.mat']), 'infoScorebootAll','sessionsSelected')



