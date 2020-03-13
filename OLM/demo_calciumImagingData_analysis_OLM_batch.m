%%%%%%%%% Calcium imaging data analysis (April-2019)-(updated July-2019) %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% README start %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This a general pipeline for performing Calcium imaging data analysis of OLM experiemnt

% Functionality of this code:
%  (1) perform basic analysis of calcium imaging data, including
%  (2) calculate the firing rate and information score information,
%  (3) identify place cells
%  (4) calcium discrimination index analysis
%  (5) compare the firing rate associated with objects
%  (6) comparison of spatial information score and firing rate
%  (7) place field analysis

%%%%%%%%%%%%%%%%%%%% README end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('/Users/suoqinjin/Documents/Yanjun_nn_revision_exp/miniscopeDataAnalysisCode_SJ_Jan2019'))

mice_saline = {'M3421F','M3423F','M3427F','M3411','M3414'};
mice_CNO = {'M3425F','M3426F','M3422F','M3424F','M3415','M3412'};
mice_rev_saline = {'M3425F_rev','M3413_rev','M3424F_rev','M3415_rev','M3412_rev'};
mice_rev_CNO = {'M3421F_rev','M3423F_rev','M3422F_rev','M3411_rev','M3414_rev'};

warning('off','all')
mice = [mice_saline,mice_CNO];
mice = [mice_rev_saline,mice_rev_CNO];

for i = 1:length(mice)
    cd('/Users/suoqinjin/Documents/Yanjun_nn_revision_exp')
    filefolder = fullfile(pwd,mice{i})
    demo_calciumImagingData_analysis_OLM(filefolder,mice{i})
    %demo_calciumImagingData_analysis_OLM_rev(filefolder,mice{i}) %
end

mice = [mice_saline,mice_CNO];
mice = [mice_rev_saline,mice_rev_CNO];
for i = 1:length(mice)
    cd('/Users/suoqinjin/Documents/Yanjun_nn_revision_exp')
    filefolder = fullfile(pwd,mice{i})
    demo_calciumImagingData_analysis_rev_placeCells(filefolder,mice{i})
    %demo_calciumImagingData_analysis_placeCells(filefolder,mice{i})
end

%% determine place cells based on the shuffling info score (info/sec or info/spike)
for i = 1:length(mouse)
    cd('/Users/suoqinjin/Documents/Yanjun_nn_revision_exp')
    filefolder = fullfile(pwd,mouse{i})
    load(fullfile(filefolder,[mouse{i},'_results_bin1_5_3SD.mat']),'infoScorebootAll','sessionCtrl')
    place_cellsAll = cell(1,size(infoScorebootAll,2));place_cellsThreshAll = [];
    
    sessionCtrl = size(infoScorebootAll,2) - 1;
    infoScoreboot = infoScorebootAll{1,sessionCtrl}; % identify place cells using info/sec
    % infoScoreboot = infoScorebootAll{2,sessionCtrl}; % identify place cells using info/spike
    
    infoScorenull = infoScoreboot(:,1); infoScoreboot = infoScoreboot(:,2:end);
    [place_cells,infoScoreThresh] = determinePlaceCells(infoScorenull,infoScoreboot);
    place_cellsAll{sessionCtrl} = place_cells;
    place_cellsThreshAll(:,sessionCtrl) = infoScoreThresh;
    save(fullfile(filefolder,['place_cells_',mouse{i},'_basedInfoSpike.mat']), 'place_cellsAll','place_cellsThreshAll')
end

%%%%%%%%%% ------- downstream analysis ------------%%%%%%%%%%
%% calcium discrimination index analysis
mice = [mice_saline,mice_CNO,mice_rev_saline,mice_rev_CNO];
for i = 1:length(mice)
    % cd('/Users/suoqinjin/Documents/NN_revision2/Yanjun_nn_revision_exp')
    filefolder = fullfile(pwd,mice{i})
    ensemble_analysis_individualMice(filefolder,mice{i})
end
mkdir results
measureNames = {'sumFiringRate','sumCountTime','sumCount','sumAmplitude','sumFiringRate_individual','sumCountTime_individual','sumCount_individual','sumAmplitude_individual'};
measureNames = {'sumFiringRate','sumCountTime','sumCount','sumAmplitude'};
% collect all the saline mice
mice = [mice_saline,mice_rev_saline];
DIall = cell2struct(cell(length(measureNames),1),measureNames)
for i = 1:length(mice)
    filefolder = fullfile(pwd,mice{i})
    load(fullfile(filefolder,'results',['discriminationIndex_',mice{i},'_radius3_allCells_sum.mat']))
    %measureNames = fieldnames(DI);
    for j = 1:length(measureNames)
        if i <= length(mice_saline)
            DIall.(char(measureNames(j))) = [DIall.(char(measureNames(j)));DI.(char(measureNames(j)))(2:end)];
        else
            DIall.(char(measureNames(j))) = [DIall.(char(measureNames(j)));DI.(char(measureNames(j)))];
        end
    end
end
DIall_saline = DIall;
save(fullfile(pwd,'results','discriminationIndex_saline.mat'), 'DIall_saline')

% collect all the CNO mice
mice = [mice_CNO,mice_rev_CNO];
DIall = cell2struct(cell(length(measureNames),1),measureNames)
for i = 1:length(mice)
    filefolder = fullfile(pwd,mice{i})
    load(fullfile(filefolder,'results',['discriminationIndex_',mice{i},'_radius3_allCells_sum.mat']))
    %  measureNames = fieldnames(DI);
    for j = 1:length(measureNames)
        if i <= length(mice_CNO)
            DIall.(char(measureNames(j))) = [DIall.(char(measureNames(j)));DI.(char(measureNames(j)))(2:end)];
        else
            DIall.(char(measureNames(j))) = [DIall.(char(measureNames(j)));DI.(char(measureNames(j)))];
        end
    end
end
DIall_CNO = DIall;
save(fullfile(pwd,'results','discriminationIndex_CNO.mat'), 'DIall_CNO')

% box plot showing the comparison of discrimination index between different groups
sessionID = 2; sessionName = 'training';
discriminationIndexComparison(DIall_saline,DIall_CNO,sessionID,sessionName)
sessionID = 3; sessionName = 'testing';
discriminationIndexComparison(DIall_saline,DIall_CNO,sessionID,sessionName)

% generate a curve showing DI with different radius
r = [2 3 4 5];
DIall_mu = zeros(length(r),3); DIall_sd = zeros(length(r),3);
for i = 1:length(r)
    filefolder = ['results_bin15_3SD_radius',num2str(r(i)),'_all_cell'];
    load(fullfile(filefolder,'discriminationIndex_CNO.mat'));
    DIall = DIall_CNO;
    DIall_mu(i,:) = mean(DIall.sumFiringRate*100);
    DIall_sd(i,:) = std(DIall.sumFiringRate*100);
end
DIall_mu_CNO = DIall_mu; DIall_sd_CNO = DIall_sd;
DIall_mu_saline = DIall_mu; DIall_sd_saline = DIall_sd;
sessionID = 2; sessionName = 'training';
discriminationIndexComparison_differentRadius(r, DIall_mu_saline,DIall_sd_saline, DIall_mu_CNO, DIall_sd_CNO,sessionID,sessionName)
sessionID = 3; sessionName = 'testing';
discriminationIndexComparison_differentRadius(r, DIall_mu_saline,DIall_sd_saline, DIall_mu_CNO, DIall_sd_CNO,sessionID,sessionName)



%% compare the firing rate associated with objects
mice = [mice_saline,mice_CNO,mice_rev_saline,mice_rev_CNO];
mice = [mice_saline,mice_rev_saline];
sumFiringRateObjectAllAll_saline = [];
for i = 1:length(mice)
    filefolder = fullfile(pwd,mice{i})
    load(fullfile(filefolder,'results',['ensembleMeasures_',mice{i},'_radius5_allCells_sum.mat']))
    sumFiringRateObjectAllAll_saline = [sumFiringRateObjectAllAll_saline;sumFiringRateObjectAll(end-2:end)];
end
mice = [mice_CNO,mice_rev_CNO];
sumFiringRateObjectAllAll_CNO = [];
for i = 1:length(mice)
    filefolder = fullfile(pwd,mice{i})
    load(fullfile(filefolder,'results',['ensembleMeasures_',mice{i},'_radius5_allCells_sum.mat']))
    sumFiringRateObjectAllAll_CNO = [sumFiringRateObjectAllAll_CNO;sumFiringRateObjectAll(end-2:end)];
end

% box plot showing the comparison of discrimination index between different groups
firingRateComparison_object(sumFiringRateObjectAllAll_saline,sumFiringRateObjectAllAll_CNO)

%% comparison of spatial information score and firing rate
sessions = {'baseline2','training','testing'};
condition = 'CNO';
if strcmpi(condition,'Saline')
    mice = [mice_saline,mice_rev_saline];
elseif strcmpi(condition,'CNO')
    mice = [mice_CNO,mice_rev_CNO];
end
InfoScoreAllAll = [];miceIDAllAll = [];MeanFiringRateAllAll = [];
for i = 1:length(mice)
    filefolder = fullfile(pwd,mice{i})
    load(fullfile(filefolder,'results',['InfoScore_',mice{i},'.mat']))
    load(fullfile(filefolder,'results',['FiringRate_',mice{i},'.mat']))
    load(fullfile(filefolder,'results',['place_cells_',mice{i},'.mat']))
    
    % segments = place_cellsAll{length(place_cellsAll)-2};
    segments = 1:size(FiringRateAll,1);
    miceIDAllAll = [miceIDAllAll; ones(length(segments),1)*i];
    
    InfoScoreAllAll = [InfoScoreAllAll;InfoScoreAll(segments,length(place_cellsAll)-1:end)];
    MeanFiringRateAllAll = [MeanFiringRateAllAll;MeanFiringRateAll(segments,length(place_cellsAll)-1:end)];
end
save(fullfile(pwd,'results',['InfoScoreAllAll_',condition,'.mat']), 'InfoScoreAllAll', 'MeanFiringRateAllAll', 'miceIDAllAll')

colors_session = []; show_pvalue = true;
SpatialInfoComparison(InfoScoreAllAll,colors_session,sessions,show_pvalue,'SI',condition) % box plot showing the distribution of info score
SpatialInfoComparison(MeanFiringRateAllAll,colors_session,sessions,show_pvalue,'event rate',condition) % box plot showing the distribution of event rate

load('InfoScoreAllAll_CNO.mat')
InfoScoreAllAll_CNO = InfoScoreAllAll;
MeanFiringRateAllAll_CNO = MeanFiringRateAllAll;
load('InfoScoreAllAll_Saline.mat')
InfoScoreAllAll_Saline = InfoScoreAllAll;
MeanFiringRateAllAll_Saline = MeanFiringRateAllAll;
conditions = {'Saline','CNO'};show_pvalue = 1; measure = 'SI';
colors_conditions = [228,26,28;55,126,184]/255;
SpatialInfoComparison_conditions(InfoScoreAllAll_Saline,InfoScoreAllAll_CNO,colors_conditions,conditions,show_pvalue,measure)

SpatialInfoComparison_conditions(MeanFiringRateAllAll_Saline,MeanFiringRateAllAll_CNO,colors_conditions,conditions,show_pvalue,'event rate')

%% compare the behavior exploration time
% behavior exploration time (10 mins)
DI_behav_CNO = [14.08 -5.24 5.08 9.19 -4.06 -4.48 12.69 5.32 -9.40-2.02 -0.73];
DI_behav_saline = [ 12.21 7.86 19.90 27.54 24.09 15.82 18.00 35.53 27.00 24.54];
sessionID = 3; sessionName = 'testing';
discriminationIndexComparison(DIall_saline,DIall_CNO,sessionID,sessionName)

% behavior exploration time (5 mins)
DI_behav_saline = [10.12 16.21 25.94 40.04 25.93 21.47 2.25 35.88 27.59 22.05];
DI_behav_CNO  = [10.00 -3.44 21.87 14.43 -10.39 -4.68 14.13 1.69 -15.53 7.40 0.51];
sessionID = 3; sessionName = 'testing';
discriminationIndexComparison(DIall_saline,DIall_CNO,sessionID,sessionName)


%% place field analysis
mice = [mice_saline,mice_CNO,mice_rev_saline,mice_rev_CNO];
for i = 1:length(mice)
    filefolder = fullfile(pwd,mice{i})
    place_field_analysis_OLM(filefolder,mice{i})
end

% collect all the place cell information
sessions = {'training','testing'};
condition = 'Saline';
if strcmpi(condition,'Saline')
    mice = [mice_saline,mice_rev_saline];
elseif strcmpi(condition,'CNO')
    mice = [mice_CNO,mice_rev_CNO];
end
% use the place cells in training
placeCellsAllAll = []; InfoPerSecAllAll = [];InfoPerSpikeAllAll = []; miceIDAllAll = [];MeanFiringRateAllAll = [];
PFresponseAllAll = []; peakERresponseAllAll = [];peakRateAllAll = [];
for i = 1:length(mice)
    filefolder = fullfile(pwd,mice{i})
    load(fullfile(filefolder,[mice{i},'_results_bin15_3SD.mat']), 'MeanFiringRateAll')
    load(fullfile(filefolder,'results',['place_cells_',mice{i},'_individualShuffle.mat']),'place_cellsAll','infoScorebootAll')
    
    segments = place_cellsAll{length(place_cellsAll)-1};
    
    placeCellsAllAll = [placeCellsAllAll; segments];
    miceIDAllAll = [miceIDAllAll; ones(length(segments),1)*i];
    
    %     InfoScoreAllAll = [InfoScoreAllAll;InfoScoreAll(segments,length(place_cellsAll):end)];
    InfoPerSecAll = [infoScorebootAll{1,length(place_cellsAll)-1}(segments,1), infoScorebootAll{1,length(place_cellsAll)}(segments,1)];
    InfoPerSecAllAll = [InfoPerSecAllAll; InfoPerSecAll];
    InfoPerSpikeAll = [infoScorebootAll{2,length(place_cellsAll)-1}(segments,1), infoScorebootAll{2,length(place_cellsAll)}(segments,1)];
    InfoPerSpikeAllAll = [InfoPerSpikeAllAll; InfoPerSpikeAll];
    MeanFiringRateAllAll = [MeanFiringRateAllAll;MeanFiringRateAll(segments,length(place_cellsAll):end)];
    
    load(fullfile(filefolder,'results',['placeFieldResponse_',mice{i},'_radius3_placeCells_basedTraining.mat']), 'PFresponseAll','peakERresponseAll','peakRateAll')
    PFresponseAllAll = [PFresponseAllAll; [PFresponseAll{length(place_cellsAll)-1}, PFresponseAll{length(place_cellsAll)}]];
    peakERresponseAllAll = [peakERresponseAllAll; [peakERresponseAll{length(place_cellsAll)-1}, peakERresponseAll{length(place_cellsAll)}]];
    peakRateAllAll = [peakRateAllAll; [peakRateAll{length(place_cellsAll)-1}, peakRateAll{length(place_cellsAll)}]];
end

InfoScoreAll = array2table([placeCellsAllAll, InfoPerSecAllAll, InfoPerSpikeAllAll, MeanFiringRateAllAll, ...
    PFresponseAllAll, peakERresponseAllAll, peakRateAllAll, miceIDAllAll], ...
    'VariableNames',{'placeCellID','ISec_training','ISec_testing','ISpike_training','ISpike_testing','FR_training','FR_testing',...
    'PF_training_obj1', 'PF_training_obj2', 'PF_training_noRes', 'PF_testing_obj1', 'PF_testing_obj2', 'PF_testing_noRes',...
    'peakER_training_obj1', 'peakER_training_obj2', 'peakER_training_noRes','peakER_testing_obj1', 'peakER_testing_obj2', 'peakER_testing_noRes',...
    'peakER_training','peakER_testing','mouse'});
save placeCellsInfoScoreAllmice_OLM_placeFieldAnalysis_basedTrainingPC_saline.mat InfoScoreAll


% use the place cells in training and testing individually
session = 'training'; in_session = 0;
session = 'testing'; in_session = 1;
placeCellsAllAll = [];  InfoPerSecAllAll = [];InfoPerSpikeAllAll = []; miceIDAllAll = [];MeanFiringRateAllAll = [];
PFresponseAllAll = []; peakERresponseAllAll = [];
for i = 1:length(mice)
    filefolder = fullfile(pwd,mice{i})
    load(fullfile(filefolder,[mice{i},'_results_bin15_3SD.mat']),  'InfoScoreAll','MeanFiringRateAll')
    load(fullfile(filefolder,'results',['place_cells_',mice{i},'_individualShuffle.mat']),'place_cellsAll','infoScorebootAll')
    
    segments = place_cellsAll{length(place_cellsAll)-1 + in_session};
    
    placeCellsAllAll = [placeCellsAllAll; segments];
    miceIDAllAll = [miceIDAllAll; ones(length(segments),1)*i];
    
    %InfoScoreAllAll = [InfoScoreAllAll;InfoScoreAll(segments,length(place_cellsAll)+in_session)];
    InfoPerSecAllAll = [InfoPerSecAllAll; infoScorebootAll{1,length(place_cellsAll)-1 + in_session}(segments,1)];
    InfoPerSpikeAllAll = [InfoPerSpikeAllAll; infoScorebootAll{2,length(place_cellsAll)-1 + in_session}(segments,1)];
    MeanFiringRateAllAll = [MeanFiringRateAllAll;MeanFiringRateAll(segments,length(place_cellsAll)+in_session)];
    
    load(fullfile(filefolder,'results',['placeFieldResponse_',mice{i},'_radius3_placeCells.mat']), 'PFresponseAll','peakERresponseAll')
    PFresponseAllAll = [PFresponseAllAll; PFresponseAll{length(place_cellsAll)-1+in_session}];
    peakERresponseAllAll = [peakERresponseAllAll; peakERresponseAll{length(place_cellsAll)-1+in_session}];
end

InfoScoreAll = array2table([placeCellsAllAll, InfoPerSecAllAll, InfoPerSpikeAllAll, MeanFiringRateAllAll, ...
    PFresponseAllAll, peakERresponseAllAll, miceIDAllAll], ...
    'VariableNames',{'placeCellID','ISec_training','ISpike_training','FR_training',...
    'PF_training_obj1', 'PF_training_obj2', 'PF_training_noRes', ...
    'peakER_training_obj1', 'peakER_training_obj2', 'peakER_training_noRes',...
    'mouse'});
save placeCellsInfoScoreAllmice_OLM_placeFieldAnalysis_individualPlaceCells_training_CNO.mat InfoScoreAll

InfoScoreAll = array2table([placeCellsAllAll, InfoPerSecAllAll, InfoPerSpikeAllAll, MeanFiringRateAllAll, ...
    PFresponseAllAll, peakERresponseAllAll, miceIDAllAll], ...
    'VariableNames',{'placeCellID','ISec_testing','ISpike_testing','FR_testing',...
    'PF_testing_obj1', 'PF_testing_obj2', 'PF_testing_noRes', ...
    'peakER_testing_obj1', 'peakER_testing_obj2', 'peakER_testing_noRes',...
    'mouse'});
save placeCellsInfoScoreAllmice_OLM_placeFieldAnalysis_individualPlaceCells_testing_CNO.mat InfoScoreAll

%% generate plots for place field analysis
% please see place_cell_analysis_OLM_basedTrainingPCs.m for details
