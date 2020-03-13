%%%%%%%%% Calcium imaging data analysis (Jan-2019)-(updated April-2019)-(updated July-2019) %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% README start %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code load neuron and behavior variable files from 'filefolder' and save results into 'filefolder' or 'filefolder/results' or 'filefolder/results/figures'  
% Please see the general pipeline for performing Calcium imaging data analysis: demo_calciumImagingData_analysis_linearTrack.m 
% for detailed description of the codes. 
% For open field, we used different parameters such as binsize, occThresh and different methods for determining the neuron firing threshold. 
%%%%%%%%%%%%%%%%%%%% README end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function demo_calciumImagingData_analysis_openField(filefolder,mice)

folderName = fullfile(filefolder,'results','figures');
if ~exist(folderName, 'dir')
    mkdir(folderName);
end
cd(filefolder)

%% open field
%% for batch analysis of multiple mice
file_neuronIndividuals = fullfile(filefolder,['neuronIndividualsf.mat'])
load(file_neuronIndividuals)
neuronIndividuals = neuronIndividualsf;
load(fullfile(filefolder,['behavIndividualsf.mat']))
behavIndividuals = behavIndividualsf;
sessionsSelected = [1 2 3];
sessionCtrl = 1;
sessions = strcat('session',cellstr(num2str([1:3]')));

% neuronfile = fullfile(filefolder,[mice,'_Neuron_filter.mat'])
% load(neuronfile)
threshfile = fullfile(filefolder,'thresh.mat')
load(threshfile)

%% for individual analysis of one mice
%% load the combined neuron data file
% addpath(genpath('/Users/suoqinjin/Documents/Yanjun_nn_revision_exp/miniscopeDataAnalysisCode_SJ_Jan2019'))
% filefolder = '/Volumes/Seagate Backup Plus Drive/NN_revision/Miniscope_OLM_preprocess/M3411';
% clear
% filefolder = pwd;
% mice = 'M3412';
% neuronfile = fullfile(filefolder,['variables_',mice,'_withTrails.mat'])
% load(neuronfile)
% file_neuronIndividuals = neuronfile
% %sessions = {'baseline1','baseline2','training','testing'};
% sessions = strcat('session',cellstr(num2str([1:length(neuronIndividuals)]')));
% load(fullfile(filefolder,['sessionsSlected_',mice,'_R1.mat']))
% sessionCtrl = sessionsSlected(1);

% neuronfile = fullfile(filefolder,[mice,'_filter_Neuron.mat'])
% load(neuronfile)
% file_neuronIndividuals = fullfile(filefolder,['neuronIndividualsf.mat'])
% neuronIndividuals = neuronIndividualsf;
% load(fullfile(filefolder,['behavIndividualsf.mat']))
% behavIndividuals = behavIndividualsf;
% sessionCtrl = 1; 
% sessions = strcat('session',cellstr(num2str([1:length(neuronIndividuals)]')));

%%
%{
%% preprocessing
% add frame infomation of neuron data if not exist
% size(neuron.C,2)
% neuron.num2read = [35981        8995        8995        8996        8995]; % make sure num2read is a vector with length(sessions)+1; first entry is total number of frames, then the number of frames in secssion 2, etc.
% isequal(neuron.num2read(1),sum(neuron.num2read(2:end)))
% save(neuronfile, 'neuron') % save neuron variable with frame information

% % split the neuron data into individual sessions
msCamId = ones(1,length(sessions));
timestamp = fullfile(filefolder,strcat('timestamp_',sessions,'.dat'));
neuronIndividuals = splittingMultiConditionNeuronData(neuron,neuronfile,msCamId,timestamp);
file_neuronIndividuals = fullfile(filefolder,['neuronIndividuals_',mice,'.mat'])
save(file_neuronIndividuals, 'neuronIndividuals','-v7.3')

% file_neuronIndividuals = fullfile(filefolder,['neuronIndividuals_new.mat'])
% load(file_neuronIndividuals)
% neuronIndividuals = neuronIndividuals_new;
% save(file_neuronIndividuals,'neuronIndividuals')
% load behav data
behavIndividuals = cell(1,length(neuronIndividuals));
for i = 1:length(neuronIndividuals)
    load(fullfile(filefolder,[sessions{i},'_Behav.mat']));
    behav0 = behav;
    t = find(diff(behav0.time)<=0);
    while ~isempty(t)
        behav0.time(t+1) = behav0.time(t)+1;
        t = find(diff(behav0.time)<=0);
    end
    %     dx = [0; diff(behav0.position(:,1))];
    %     dy = [0; diff(behav0.position(:,2))];
    %     behav0.speed = sqrt((dx).^2+(dy).^2)/behav0.dt;
    %     behav0.speed = smoothts(behav0.speed','b',ceil(1/behav0.dt/5));
    behavIndividuals{i} = behav0;
end
file_behavIndividuals = fullfile(filefolder,['behavIndividuals_',mice,'.mat'])
save(file_behavIndividuals, 'behavIndividuals')

% % % % delete the first data point in the time information if the first frame is filtered out
% % % for i = 1:length(neuronIndividuals)
% % %     neuronIndividuals{i}.time(1) = [];
% % %     behavIndividuals{i}.time(1) = [];
% % % end

% visualization
behavPlot_parallel(behavIndividuals,sessions)
% plottingTrace(neuron,1:size(neuron.C,1),1)
%}

%% ------------ analysis start: ------------------
%% determine the neuron firing threshold
temp = 'S';
% thresh = determiningFiringEventThresh(neuron,temp,'sd');

%% calculate the firing rate and information score information
%% incorporate the firing rate information and then plot the final spatial coding heat map
% binsize = 10; countTimeThresh = 0.2;
binsize = 2; countTimeThresh = 0.1;% open filed
% [FiringRateAll,CountAll,~,CountTimeAll1,MeanFiringRateAll,InfoScoreAll,InfoSpikeAll] = calculating_FR_IS_parallel(neuronIndividuals,behavIndividuals,binsize,temp,thresh,countTimeThresh);
% save(fullfile(filefolder,'results',['FiringRate_',mice,'.mat']), 'FiringRateAll','CountAll','CountTimeAll1','MeanFiringRateAll')
% save(fullfile(filefolder,'results',['InfoScore_',mice,'.mat']), 'InfoScoreAll','InfoSpikeAll')

% comparison of information score and mean firing rate of all cells

% colors_session = []; show_pvalue = true;
% SpatialInfoComparison(InfoScoreAll(:,2:end),colors_session,sessions,show_pvalue,'SI')
% SpatialInfoComparison(MeanFiringRateAll(:,2:end),colors_session,sessions,show_pvalue,'event rate')
% segments = 1:size(neuronIndividuals{1}.C,1);
% threshSpatial = determiningSpatialThresh(FiringRateAll,segments)
% threshSpatial = 1;
% overlayFiringRate(neuronIndividuals,behavIndividuals,FiringRateAll,segments(1:10),thresh,temp,threshSpatial,strcat(mice,'_',sessions,'_all'))
%
%close all
%% identify place cells
occThresh = 0.1; nboot = 1000;
% randomly generate the dealt t for perpute the spike
deltaTall = randi([10,890],nboot,1)*1000; % unit: ms
% neuron_lowFR = find(sum(MeanFiringRateAll(:,2:end) < 0.1,2) >= 2);
neuron_lowFR = [];
experiment = 'openField';
[place_cellsAll,place_cellsThreshAll,infoScorebootAll] = permutingSpike_parallel(file_neuronIndividuals,behavIndividuals,thresh,temp,deltaTall,occThresh,nboot,binsize,sessions,sessionCtrl,neuron_lowFR,experiment);
% [place_cells,TinfoPerSecond] = permutingSpike(file_neuronIndividuals,sessionIndex,neuron0,behav0,thresh,'S',deltaTall,occThresh,nboot,binsize);
save(fullfile(filefolder,'results',['place_cells_',mice,'.mat']), 'place_cellsAll','place_cellsThreshAll','infoScorebootAll','sessionCtrl')


% comparison of information score and mean firing rate of all cells

% colors_session = []; show_pvalue = true;
% SpatialInfoComparison(InfoScoreAll(:,2:end),colors_session,sessions,show_pvalue,'SI')
% SpatialInfoComparison(MeanFiringRateAll(:,2:end),colors_session,sessions,show_pvalue,'event rate')
% segments = 1:size(neuronIndividuals{1}.C,1);
% threshSpatial = determiningSpatialThresh(FiringRateAll,segments)
% threshSpatial = 1;
% overlayFiringRate(neuronIndividuals,behavIndividuals,FiringRateAll,segments(1:10),thresh,temp,threshSpatial,strcat(mice,'_',sessions,'_all'))
%

% clear neuron neuronIndividuals behavIndividuals
% save results_measures_M3411_bin10_10peak.mat

%{
% comparison of information score and mean firing rate for place cells
segments = place_cellsAll{sessionCtrl};

colors_session = []; show_pvalue = true;
SpatialInfoComparison(InfoScoreAll(segments,2:end),colors_session,sessions,show_pvalue,'SI') % box plot showing the distribution of info score
SpatialInfoComparison(MeanFiringRateAll(segments,2:end),colors_session,sessions,show_pvalue,'event rate') % box plot showing the distribution of event rate
SpatialInfoComparison(InfoSpikeAll(segments,2:end),colors_session,sessions,show_pvalue,'SIspike') 

%threshSpatial = determiningSpatialThresh(FiringRateAll,segments)
threshSpatial = 3;
overlayFiringRate(neuronIndividuals,behavIndividuals,FiringRateAll,segments(1:min(20,length(segments))),thresh,temp,threshSpatial,strcat(mice,'_',sessions,'_placeCells'))
close all
clear neuron neuronIndividuals behavIndividuals behav behav0

save(fullfile(filefolder,[mice,'_results_bin10_10peak_basedInfoSec.mat']))
% save(fullfile(filefolder,[mice,'_results_bin10_10peak.mat'])) % linear track
% save(fullfile(filefolder,[mice,'_results_bin15_3SD.mat'])) % OLM
%save(fullfile(filefolder,[mice,'_results_bin1_5_10peak.mat'])) % open field
%save(fullfile(filefolder,[mice,'_results_bin1_5_3SD.mat'])) % open field
%}





