function ensemble_analysis_individualMice(filefolder,mice)
%%%%%%%%%%%%%%%%%%%% README start %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code calculates the ensemble information associated with objects in OLM experiemnt, including firing rate, exploration time and number of events
% You might need to change the parameter 'radius' to calculate the ensemble information associated with the object.
%%%%%%%%%%%%%%%%%%%% README end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

resultsfile = fullfile(filefolder,[mice,'_results_bin15_3SD.mat'])
load(resultsfile)
% load(fullfile(filefolder,'results',['FiringRate_',mice,'.mat']))
% load(fullfile(filefolder,'results',['place_cells_',mice,'.mat']))
load(fullfile(filefolder,['behavIndividuals_',mice,'.mat']))

% folderName = fullfile(filefolder,'results')
% if ~exist(folderName, 'dir')
%     mkdir(folderName);
% end

%% --------start to calculate ensemble information ---------------------------
% behav = behavIndividuals{4};firingRate = FiringRateAll(:,4);binsize = 15; pos = [(behav.object(:,1)-behav.ROI(1))*behav.trackLength/behav.ROI(3),behav.object(:,2)*behav.trackLength/behav.ROI(3)];
% measure = 'event rate';colorscale = [0 0.1];plotting = 1;addindex = 0;objname = {'Orig','New'};
% [sumFiringRate,sumFiringRate_conv,sumFiringRate_smooth,sumFiringRateObject,posObjects,sumFiringRateObject_individual,radius] = comparingFiringRateSingleConditionMultiObjects(firingRate,binsize,CountTime, pos,measure,colorscale,plotting,addindex,objname);

% segments = place_cellsAll{length(place_cellsAll)-1}; % only use place cells in baseline 2 for analysis
segments = [];
plotting = 1;addindex = 0;objname = {'Object1','Object2'};smooth_index = 1;

radius = 3; %individual neuron's event rate associated with the object was quantified as the average event rate across bins within (radius-1)bins surrounding the object.

measures = 'event rate';colorscale = [0 50];
[sumFiringRateAll,sumFiringRate_convAll,sumFiringRate_smoothingAll,sumFiringRateObjectAll,sumFiringRateObject_individualAll,posObjectsAll] = ...
    comparingFiringRateMultiConditionMultiObjects(behavIndividuals,FiringRateAll,segments,CountTimeAll,binsize,binInfoAll,radius,measures,plotting,smooth_index,sessions,colorscale,addindex,objname);

measures = 'occupation time';colorscale = [0 2];
[sumCountTimeAll,sumCountTime_convAll,sumCountTime_smoothingAll,sumCountTimeObjectAll,sumCountTimeObject_individualAll,posObjectsAll] = ...
    comparingFiringRateMultiConditionMultiObjects(behavIndividuals,CountTimeAll,segments,CountTimeAll,binsize,binInfoAll,radius,measures,plotting,smooth_index,sessions,colorscale,addindex,objname);

measures = 'count';colorscale = [0 0.4];
[sumCountAll,sumCount_convAll,sumCount_smoothingAll,sumCountObjectAll,sumCountObject_individualAll,posObjectsAll] = ...
    comparingFiringRateMultiConditionMultiObjects(behavIndividuals,CountAll,segments,CountTimeAll,binsize,binInfoAll,radius,measures,plotting,smooth_index,sessions,colorscale,addindex,objname);

measures = 'amplitude';colorscale = [0 0.3]; % NB: the calculated amplitude might be not correct
[sumAmplitudeAll,sumAmplitude_convAll,sumAmplitude_smoothingAll,sumAmplitudeObjectAll,sumAmplitudeObject_individualAll,posObjectsAll] = ...
    comparingFiringRateMultiConditionMultiObjects(behavIndividuals,AmplitudeAll,segments,CountTimeAll,binsize,binInfoAll,radius,measures,plotting,smooth_index,sessions,colorscale,addindex,objname);

save(fullfile(filefolder,'results',['ensembleMeasures_',mice,'_radius',num2str(radius),'_allCells_sum.mat']), 'sumFiringRateAll','sumFiringRate_convAll','sumFiringRate_smoothingAll','sumFiringRateObjectAll','sumFiringRateObject_individualAll',...
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

save(fullfile(filefolder,'results',['discriminationIndex_',mice,'_radius',num2str(radius),'_allCells_sum.mat']), 'DI')


