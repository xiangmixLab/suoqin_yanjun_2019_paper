%%%%%%%%% Calcium imaging data analysis (Jan-2019)-(updated April-2019)-(updated July-2019) %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% README start %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is a general pipeline for performing jackknife or bootstrap sampling analysis, including linear track and open field
% This code load neuron and behavior variable files from 'filefolder' and save results into 'filefolder'
% Functionality of this code: recalculate info score using jackknife or bootstrap sampling method
%%%%%%%%%%%%%%%%%%%% README end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function demo_calciumImagingData_analysis_linearTrack_subseting(filefolder,mice, method_sampling)
if ~exist('method_sampling','var') || isempty(method_sampling)
    method_sampling = 'jackknife';
end

folderName = fullfile(filefolder,'results','figures');
if ~exist(folderName, 'dir')
    mkdir(folderName);
end
cd(filefolder)
%% CNO-linear track
% load combined neuron variable for determining the firing threshold
% load individual neuron file and behavior file
neuronfile = fullfile(filefolder,['variables_',mice,'_withTrails.mat'])
load(neuronfile)
file_neuronIndividuals = neuronfile
sessions = strcat('session',cellstr(num2str([1:length(neuronIndividuals)]')));
load(fullfile(filefolder,['sessionsSlected_',mice,'_R1.mat']))
sessionUsed = sessionsSlected;
%% saline control-linear track
% neuronfile = fullfile(filefolder,[mice,'_filter_Neuron.mat'])
% load(neuronfile)
% file_neuronIndividuals = fullfile(filefolder,['neuronIndividualsf.mat'])
% load(file_neuronIndividuals)
% neuronIndividuals = neuronIndividualsf;
% load(fullfile(filefolder,['behavIndividualsf.mat']))
% behavIndividuals = behavIndividualsf;
% sessionCtrl = 1;
% sessions = strcat('session',cellstr(num2str([1:length(neuronIndividuals)]')));
% sessionsSlected = [1 2 3];


%% ------------ analysis start: ------------------
%% determine the neuron firing threshold
temp = 'S';
thresh = determiningFiringEventThresh(neuron,temp);

%% calculate the firing rate and information score information
%% incorporate the firing rate information and then plot the final spatial coding heat map
binsize = 10; countTimeThresh = 0.2;
% %binsize = 1.5; countTimeThresh = 0.1;% open filed
% [FiringRateAll,CountAll,~,CountTimeAll1,MeanFiringRateAll,InfoScoreAll,InfoSpikeAll] = calculating_FR_IS_parallel(neuronIndividuals,behavIndividuals,binsize,temp,thresh,countTimeThresh);
% save(fullfile(filefolder,'results',['FiringRate_',mice,'.mat']), 'FiringRateAll','CountAll','CountTimeAll1','MeanFiringRateAll')
% save(fullfile(filefolder,'results',['InfoScore_',mice,'.mat']), 'InfoScoreAll','InfoSpikeAll')
%close all
%% identify place cells
% occThresh = 0.2; nboot = 1000;
% % randomly generate the dealt t for perpute the spike
% deltaTall = randi([10,890],nboot,1)*1000; % unit: ms
% % neuron_lowFR = find(sum(MeanFiringRateAll(:,2:end) < 0.1,2) >= 2);
% neuron_lowFR = [];
%
% [place_cellsAll,place_cellsThreshAll,infoScorebootAll] = permutingSpike_parallel(file_neuronIndividuals,behavIndividuals,thresh,temp,deltaTall,occThresh,nboot,binsize,sessions,sessionCtrl,neuron_lowFR);
% % [place_cells,TinfoPerSecond] = permutingSpike(file_neuronIndividuals,sessionIndex,neuron0,behav0,thresh,'S',deltaTall,occThresh,nboot,binsize);
% save(fullfile(filefolder,'results',['place_cells_',mice,'.mat']), 'place_cellsAll','place_cellsThreshAll','infoScorebootAll')
experiment = 'linearTrack';
occThresh = 0.2; nboot = 1000;
nbin = 10; % each track-running session was divided into 10 equal-duration sub-sessions
if strcmp(mice,'1stCA1')
    infoScorebootAll = subsetingSpike_parallel_nooverlap2(file_neuronIndividuals,behavIndividuals,thresh,temp,occThresh,nboot,binsize,sessionUsed,nbin, experiment, method_sampling);
else
    infoScorebootAll = subsetingSpike_parallel_nooverlap(file_neuronIndividuals,behavIndividuals,thresh,temp,occThresh,nboot,binsize,sessionUsed,nbin, experiment, method_sampling);
end
save(fullfile(filefolder,['infoScorebootAll_subsetingBiningTime_',num2str(nboot),mice,'.mat']), 'infoScorebootAll','sessionsSlected')





