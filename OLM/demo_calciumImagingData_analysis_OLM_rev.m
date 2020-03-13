%%%%%%%%% Calcium imaging data analysis (Jan-2019)-(updated July-2019) %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% README start %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code load neuron and behavior variable files from 'filefolder' and save results into 'filefolder' or 'filefolder/results' or 'filefolder/results/figures'
% rev experiments only have three sessions compared to the four sessions in other OLM experiments
%%%%%%%%%%%%%%%%%%%% README end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function demo_calciumImagingData_analysis_OLM_rev(filefolder,mice)

folderName = fullfile(filefolder,'results','figures');
if ~exist(folderName, 'dir')
    mkdir(folderName);
end
%% load the combined neuron data file
cd(filefolder)
neuronfile = fullfile(filefolder,['further_processed_neuron_extraction_final_results.mat'])
load(neuronfile)
sessions = {'baseline2','training','testing'};

for j = 1:size(neuron.C,1)
    [pks,locs] = findpeaks(neuron.C(j,:));
    neuron.S(j,:) = zeros(size(neuron.C(j,:)));
    neuron.S(j,locs) = pks;
end
neuronfile = fullfile(filefolder,[mice,'_neuron_finalResults.mat'])
save(neuronfile,'neuron')
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

%% ------------ analysis start: ------------------
%% determine the neuron firing threshold
temp = 'S'; method = 'sd';
thresh = determiningFiringEventThresh(neuron,temp,method);

%% calculate the firing rate and information score information
%% incorporate the firing rate information and then plot the final spatial coding heat map
binsize = 15; % the binsize is dependent on the tracklength used in the extracted behavior data
countTimeThresh = 0.1; % set firing rate to be zero if the occupation time is less than countTimeThresh
experiment = 'OLM'; % linear track, open field and OLM will use different functions to calculate firing rate maps
[FiringRateAll,CountAll,CountTimeAll,CountTimeAll1,MeanFiringRateAll,InfoScoreAll,InfoSpikeAll,AmplitudeAll,binInfoAll] = calculating_FR_IS_parallel(neuronIndividuals,behavIndividuals,binsize,temp,thresh,countTimeThresh, experiment);
save(fullfile(filefolder,'results',['FiringRate_',mice,'.mat']), 'FiringRateAll','CountAll','CountTimeAll','CountTimeAll1','MeanFiringRateAll','AmplitudeAll','binInfoAll')
save(fullfile(filefolder,'results',['InfoScore_',mice,'.mat']), 'InfoScoreAll')

% comparison of information score and mean firing rate of all cells

% colors_session = []; show_pvalue = true;
% SpatialInfoComparison(InfoScoreAll(:,2:end),colors_session,sessions,show_pvalue,'SI')
% SpatialInfoComparison(MeanFiringRateAll(:,2:end),colors_session,sessions,show_pvalue,'event rate')
% segments = 1:size(neuronIndividuals{1}.C,1);
% threshSpatial = determiningSpatialThresh(FiringRateAll,segments)
% threshSpatial = 1;
% overlayFiringRate(neuronIndividuals,behavIndividuals,FiringRateAll,segments(1:10),thresh,temp,threshSpatial,strcat(mice,'_',sessions,'_all'))
%
close all
%% identify place cells
occThresh = 0.1; nboot = 1000;
% randomly generate the dealt t for perpute the spike
deltaTall = randi([10,590],nboot,1)*1000; % unit: ms; in total 600s
neuron_lowFR = [];
experiment = 'OLM';
[place_cellsAll,place_cellsThreshAll,infoScorebootAll] = permutingSpike_parallel(file_neuronIndividuals,behavIndividuals,thresh,temp,deltaTall,occThresh,nboot,binsize,sessions,sessionCtrl,neuron_lowFR,experiment);
save(fullfile(filefolder,'results',['place_cells_',mice,'.mat']), 'place_cellsAll','place_cellsThreshAll','infoScorebootAll')

% comparison of information score and mean firing rate for place cells
segments = place_cellsAll{1};

colors_session = []; show_pvalue = true;
SpatialInfoComparison(InfoScoreAll(segments,2:end),colors_session,sessions,show_pvalue,'SI') % box plot showing the distribution of info score
SpatialInfoComparison(MeanFiringRateAll(segments,2:end),colors_session,sessions,show_pvalue,'event rate') % box plot showing the distribution of event rate

threshSpatial = determiningSpatialThresh(FiringRateAll,segments)
threshSpatial = 1;
overlayFiringRate(neuronIndividuals,behavIndividuals,FiringRateAll,segments,thresh,temp,threshSpatial,strcat(mice,'_',sessions,'_placeCells'))

close all
clear neuron neuronIndividuals behavIndividuals behav behav0
% save results_measures_M3411_bin10_10peak.mat
% save(fullfile(filefolder,[mice,'_results_bin10_10peak.mat']))
save(fullfile(filefolder,[mice,'_results_bin15_3SD.mat']))



