%%%%%%%%% Calcium imaging data analysis (Jan-2019) %%%%%%%%%%%%%
function demo_calciumImagingData_analysis_rev_placeCells(filefolder,mice)
%% load the combined neuron data file
% addpath(genpath('/Users/suoqinjin/Documents/Yanjun_nn_revision_exp/miniscopeDataAnalysisCode_SJ_Jan2019'))
% filefolder = '/Volumes/Seagate Backup Plus Drive/NN_revision/Miniscope_OLM_preprocess/M3411';
% clear
% filefolder = pwd;
% mice = 'M3412';
folderName = fullfile(filefolder,'results','figures');
if ~exist(folderName, 'dir')
    mkdir(folderName);
end
cd(filefolder)
%{
neuronfile = fullfile(filefolder,['further_processed_neuron_extraction_final_results.mat'])
load(neuronfile)
sessions = {'baseline2','training','testing'};

for j = 1:size(neuron.C,1)
    %         [pks,locs] = findpeaks(neuron.C(j,:),'MinPeakHeight',thresh_amp(j));
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

neuronfile = fullfile(filefolder,[mice,'_neuron_finalResults.mat'])
load(neuronfile)
file_neuronIndividuals = fullfile(filefolder,['neuronIndividuals_',mice,'.mat'])
file_behavIndividuals = fullfile(filefolder,['behavIndividuals_',mice,'.mat'])
load(file_behavIndividuals)

sessions = {'baseline2','training','testing'};
sessionsSelected = [1 2 3];
% sessionCtrl = sessionsSlected(1);

%% ------------ analysis start: ------------------
%% determine the neuron firing threshold
temp = 'S'; method = 'sd';
thresh = determiningFiringEventThresh(neuron,temp,method);

%% calculate the firing rate and information score information
%% incorporate the firing rate information and then plot the final spatial coding heat map
binsize = 15; countTimeThresh = 0.1;
%{
experiment = 'OLM';
[FiringRateAll,CountAll,CountTimeAll,CountTimeAll1,MeanFiringRateAll,InfoSecondAll,InfoSpikeAll,AmplitudeAll,binInfoAll] ...
    = calculating_FR_IS_parallel(neuronIndividuals,behavIndividuals,binsize,temp,thresh,countTimeThresh,experiment);
save(fullfile(filefolder,'results',['FiringRate_',mice,'.mat']), 'FiringRateAll','CountAll','CountTimeAll','CountTimeAll1','MeanFiringRateAll','AmplitudeAll','binInfoAll')
save(fullfile(filefolder,'results',['InfoScore_',mice,'.mat']), 'InfoSecondAll','InfoSpikeAll')
%}
% comparison of information score and mean firing rate of all cells

% colors_session = []; show_pvalue = true;
% SpatialInfoComparison(InfoScoreAll(:,2:end),colors_session,sessions,show_pvalue,'SI')
% SpatialInfoComparison(MeanFiringRateAll(:,2:end),colors_session,sessions,show_pvalue,'event rate')
% segments = 1:size(neuronIndividuals{1}.C,1);
% threshSpatial = determiningSpatialThresh(FiringRateAll,segments)
% threshSpatial = 1;
% overlayFiringRate(neuronIndividuals,behavIndividuals,FiringRateAll,segments(1:10),thresh,temp,threshSpatial,strcat(mice,'_',sessions,'_all'))
%
% close all
%% identify place cells
% occThresh = 1; nboot = 100;
occThresh = 0.1; nboot = 1000;
% randomly generate the dealt t for perpute the spike
deltaTall = randi([10,590],nboot,1)*1000; % unit: ms; in total 600s
% file_neuronIndividuals = ['neuronIndividuals_',mice,'.mat'];

% [place_cellsAll,place_cellsThreshAll] = permutingSpike_parallel(file_neuronIndividuals,behavIndividuals,thresh,temp,deltaTall,occThresh,nboot,binsize,sessions);
% % [place_cells,TinfoPerSecond] = permutingSpike(file_neuronIndividuals,sessionIndex,neuron0,behav0,thresh,'S',deltaTall,occThresh,nboot,binsize);
% save(fullfile(filefolder,'results',['place_cells_',mice,'.mat']), 'place_cellsAll','place_cellsThreshAll')
neuron_lowFR = []; experiment = 'OLM';
[place_cellsAll,place_cellsThreshAll,infoScorebootAll] = permutingSpike_parallel(file_neuronIndividuals,behavIndividuals,thresh,temp,deltaTall,occThresh,nboot,binsize,sessions,sessionsSelected,neuron_lowFR,experiment);
save(fullfile(filefolder,'results',['place_cells_',mice,'_individualShuffle.mat']), 'place_cellsAll','place_cellsThreshAll','infoScorebootAll')
% clear neuron neuronIndividuals behavIndividuals
% save results_measures_M3411_bin10_10peak.mat
%{
%% --------start to calculate ensemble information ---------------------------
% behav = behavIndividuals{4};firingRate = FiringRateAll(:,4);binsize = 15; pos = [(behav.object(:,1)-behav.ROI(1))*behav.trackLength/behav.ROI(3),behav.object(:,2)*behav.trackLength/behav.ROI(3)];
% measure = 'event rate';colorscale = [0 0.1];plotting = 1;addindex = 0;objname = {'Orig','New'};
% [sumFiringRate,sumFiringRate_conv,sumFiringRate_smooth,sumFiringRateObject,posObjects,sumFiringRateObject_individual,radius] = comparingFiringRateSingleConditionMultiObjects(firingRate,binsize,CountTime, pos,measure,colorscale,plotting,addindex,objname);
segments = place_cellsAll{1}; % only use place cells in baseline 2 for analysis
plotting = 1;addindex = 0;objname = {'Object1','Object2'};smooth_index = 1;
measures = 'event rate';colorscale = [0 0.2];
[sumFiringRateAll,sumFiringRate_convAll,sumFiringRate_smoothingAll,sumFiringRateObjectAll,sumFiringRateObject_individualAll,posObjectsAll] = ...
    comparingFiringRateMultiConditionMultiObjects(behavIndividuals,FiringRateAll,segments,CountTimeAll,binsize,binInfoAll,measures,plotting,smooth_index,sessions,colorscale,addindex,objname);

measures = 'occupation time';colorscale = [0 2];
[sumCountTimeAll,sumCountTime_convAll,sumCountTime_smoothingAll,sumCountTimeObjectAll,sumCountTimeObject_individualAll,posObjectsAll] = ...
    comparingFiringRateMultiConditionMultiObjects(behavIndividuals,CountTimeAll,segments,CountTimeAll,binsize,binInfoAll,measures,plotting,smooth_index,sessions,colorscale,addindex,objname);

measures = 'count';colorscale = [0 0.4];
[sumCountAll,sumCount_convAll,sumCount_smoothingAll,sumCountObjectAll,sumCountObject_individualAll,posObjectsAll] = ...
    comparingFiringRateMultiConditionMultiObjects(behavIndividuals,CountAll,segments,CountTimeAll,binsize,binInfoAll,measures,plotting,smooth_index,sessions,colorscale,addindex,objname);

measures = 'amplitude';colorscale = [0 0.3];
[sumAmplitudeAll,sumAmplitude_convAll,sumAmplitude_smoothingAll,sumAmplitudeObjectAll,sumAmplitudeObject_individualAll,posObjectsAll] = ...
    comparingFiringRateMultiConditionMultiObjects(behavIndividuals,CountAll,segments,CountTimeAll,binsize,binInfoAll,measures,plotting,smooth_index,sessions,colorscale,addindex,objname);

save(fullfile(filefolder,'results',['ensembleMeasures_',mice,'.mat']), 'sumFiringRateAll','sumFiringRate_convAll','sumFiringRate_smoothingAll','sumFiringRateObjectAll','sumFiringRateObject_individualAll',...
    'sumCountTimeAll','sumCountTime_convAll','sumCountTime_smoothingAll','sumCountTimeObjectAll','sumCountTimeObject_individualAll',...
    'sumCountAll','sumCount_convAll','sumCount_smoothingAll','sumCountObjectAll','sumCountObject_individualAll',...
    'sumAmplitudeAll','sumAmplitude_convAll','sumAmplitude_smoothingAll','sumAmplitudeObjectAll','sumAmplitudeObject_individualAll','posObjectsAll')
% calculate discrimination index
DI.sumFiringRate = calculate_discriminationIndex(sumFiringRateObjectAll);
DI.sumCountTime = calculate_discriminationIndex(sumCountTimeObjectAll);
DI.sumCount = calculate_discriminationIndex(sumCountObjectAll);
DI.sumAmplitude = calculate_discriminationIndex(sumAmplitudeObjectAll);

DI.sumFiringRate_individual = calculate_discriminationIndex(sumFiringRateObject_individualAll);
DI.sumCountTime_individual = calculate_discriminationIndex(sumCountTimeObject_individualAll);
DI.sumCount_individual = calculate_discriminationIndex(sumCountObject_individualAll);
DI.sumAmplitude_individual = calculate_discriminationIndex(sumAmplitudeObject_individualAll);

save(fullfile(filefolder,'results',['discriminationIndex_',mice,'.mat']), 'DI')
%}

% save(fullfile(filefolder,[mice,'_results_bin10_10peak.mat']),  'FiringRateAll','CountAll','CountTimeAll','CountTimeAll1','MeanFiringRateAll','binInfoAll',...
%    'InfoScoreAll', 'place_cellsAll','place_cellsThreshAll','DI',...
%    'sumFiringRateAll','sumFiringRate_convAll','sumFiringRate_smoothingAll','sumFiringRateObjectAll',...
%     'sumCountTimeAll','sumCountTime_convAll','sumCountTime_smoothingAll','sumCountTimeObjectAll',...
%     'sumCountAll','sumCount_convAll','sumCount_smoothingAll','sumCountObjectAll',...
%     'sumAmplitudeAll','sumAmplitude_convAll','sumAmplitude_smoothingAll','sumAmplitudeObjectAll','posObjectsAll')
%{
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
%}
% %% find the peak and then calculate the amplitude
% thresh_amp = max(neuron.C,[],2)*0.1;
% figure;findpeaks(neuron.C(5,:),'MinPeakHeight',thresh_amp(5))
% amp_avg_C = zeros(size(neuronIndividuals{1}.C,1),length(neuronIndividuals));
% for i = 1:length(neuronIndividuals)
%     neuron0 = neuronIndividuals{i};
%     for j = 1:size(neuron0.C,1)
%         [pks,locs,w,p] = findpeaks(neuron0.C(j,:),'MinPeakHeight',thresh_amp(j));
%         amp_avg_C(j,i) = mean(pks);
%     end
% end
% figure; boxplot(amp_avg_C)



