%%%%%%%%% Calcium imaging data analysis (April-2019)-(updated July-2019) %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% README start %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Please see/use the general pipeline for performing Calcium imaging data analysis, demo_calciumImagingData_analysis_linearTrack_batch.m
% The difference of this pipeline with
% demo_calciumImagingData_analysis_linearTrack_batch.m is just different function names used

% Functionality of this code:
% (1) perform basic analysis of calcium imaging data, including
%    (i) calculate the firing rate and information score information,
%    (ii) identify place cells
%    (iii) assemble all the calculated information into one table for various downstream analysis

% (2) jackknife/boostrap analysis, calculate the percent of each category by varying the threshold of pvalues from thesampling analysis
%%%%%%%%%%%%%%%%%%%% README end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('/Users/suoqinjin/Documents/Yanjun_nn_revision_exp/miniscopeDataAnalysisCode_SJ_Jan2019'))

warning('off','all')
mice_saline_open = {'M3411','M3414','M3421F','M3422F'};
mice_CNO_open = {'M3412','M3413','M3424F','M3425F'};
mouse = mice_saline_open;
for i = 1:length(mouse)
    cd('/Users/suoqinjin/Documents/NN_revision2/Yanjun_openField_results')
    filefolder = fullfile(pwd,mouse{i})
    demo_calciumImagingData_analysis_openField(filefolder,mouse{i})
end

%% determine place cells based on the shuffling info score (info/sec or info/spike)
mouse = {'M3244F','M3321','M3322','M3243','M3323','1stCA1'};
for i = 1:length(mouse)
    cd('/Users/suoqinjin/Documents/NN_revision2/Yanjun_openField_results')
    filefolder = fullfile(pwd,mouse{i})
    load(fullfile(filefolder,[mouse{i},'_results_bin1_5_3SD.mat']),'infoScorebootAll','sessionCtrl')
    place_cellsAll = cell(1,size(infoScorebootAll,2));place_cellsThreshAll = [];
    
    infoScoreboot = infoScorebootAll{1,sessionCtrl}; % identify place cells using info/sec
    % infoScoreboot = infoScorebootAll{2,sessionCtrl}; % identify place cells using info/spike
    
    infoScorenull = infoScoreboot(:,1); infoScoreboot = infoScoreboot(:,2:end);
    [place_cells,infoScoreThresh] = determinePlaceCells(infoScorenull,infoScoreboot);
    place_cellsAll{sessionCtrl} = place_cells;
    place_cellsThreshAll(:,sessionCtrl) = infoScoreThresh;
    save(fullfile(filefolder,['place_cells_',mouse{i},'_basedInfoSpike.mat']), 'place_cellsAll','place_cellsThreshAll')
end

%% assemble all the information (e.g., info score, firing rate, mouse ID and place cell ID) into a table InfoScoreAll
sessions = {'Ctrl','CNO','PCtrl'};
InfoPerSecAllAll = [];miceIDAllAll = [];MeanFiringRateAllAll = [];
placeCellsAllAll=[];InfoPerSpikeAllAll = [];
for i = 1:length(mouse)
    filefolder = fullfile(pwd,mouse{i})
    load(fullfile(filefolder,[mouse{i},'_results_bin2_3SDpeak.mat']),'InfoScoreAll','InfoSpikeAll','MeanFiringRateAll','sessionCtrl')
    load(fullfile(filefolder,['place_cells_',mouse{i},'_basedInfoSpike.mat']),'place_cellsAll')
    if ~exist('sessionsSlected','var') || isempty(sessionsSlected), sessionsSlected = [1 2 3];end
    segments = place_cellsAll{sessionCtrl};
    %segments = 1:size(FiringRateAll,1);
    placeCellsAllAll = [placeCellsAllAll; segments];
    miceIDAllAll = [miceIDAllAll; ones(length(segments),1)*i];
    InfoPerSecAllAll = [InfoPerSecAllAll;InfoScoreAll(segments,sessionsSlected+1)];
    InfoPerSpikeAllAll = [InfoPerSpikeAllAll;InfoSpikeAll(segments,sessionsSlected+1)];
    MeanFiringRateAllAll = [MeanFiringRateAllAll;MeanFiringRateAll(segments,sessionsSlected+1)];
end
InfoScoreAll = array2table([placeCellsAllAll, InfoPerSecAllAll, InfoPerSpikeAllAll, MeanFiringRateAllAll, miceIDAllAll], ...
    'VariableNames',{'placeCellID','ISec_Ctrl','ISec_CNO','ISec_PCtrl','ISpike_Ctrl','ISpike_CNO','ISpike_PCtrl',...
    'FR_Ctrl','FR_CNO','FR_PCtrl','mouse'});
save placeCellsInfoScoreAllmice_openField_CNO_basedInfoSec_binsize1_5_10peak.mat InfoScoreAll

%% calculate the SI using subset of cells (jackknife sampling)
for i = 1:length(mouse)
    cd('/Users/suoqinjin/Documents/NN_revision2/Yanjun_openField_results')
    filefolder = fullfile(pwd,mouse{i})
    method_sampling = 'jackknife';
    % recalculate info score using jackknife or bootstrap sampling method
    demo_calciumImagingData_analysis_openField_subseting(filefolder,mouse{i},method_sampling)
    % load the recalculated info score
    load(fullfile(filefolder,['infoScorebootAll_subsetingBiningTime_5_1000',mouse{i},'.mat']), 'infoScorebootAll','sessionsSlected')
    
    % calculate pvalues of each neuron based on jackknife sampling
    if strcmpi(method_sampling,'jackknife')
        [Qstatistics, Pvalues] = distributionTest_jackknife(infoScorebootAll,sessionsSelected);
        save(fullfile(filefolder,[mouse{i},'_distributionTest_jackknife_nboot10BiningTime_bitSpikeYanjun.mat']),'Qstatistics','Pvalues')
    elseif strcmpi(method_sampling,'bootstrap')
        [D12, D32, P12, P32] = distributionTest_bootstrap(infoScorebootAll,sessionsSelected);
        save(fullfile(filefolder,[mouse{i},'_distributionTest_bootstrap_nboot1000_5BiningTime_bitSpikeYanjun.mat']),'D12','D32','P12','P32')
    end
end


%% update the InfoScoreAll variable by adding jackknife sampling analysis
%%% new_individualShuffle_yanjun05015
% load the InfoScoreAll variable that contain all the informatin (e.g., info score, mice id) of all mice
load('/Users/suoqinjin/Documents/NN_revision2/Yanjun_openField_results/new_individualShuffle_yanjun0515/placeCellsInfoScoreAllmice_CNOsec_final.mat')
load('/Users/suoqinjin/Documents/NN_revision2/Yanjun_openField_results/new_individualShuffle_yanjun0515/placeCellsInfoScoreAllmice_salinesec_final.mat')
PvaluesAll = [];QstatisticsAll = [];peakFRAll = [];
for i = 1:length(mouse)
    filefolder = fullfile(pwd,mouse{i})
    load(fullfile(filefolder,[mouse{i},'_distributionTest_jackknife_nboot5BiningTime.mat']),'Qstatistics','Pvalues')
    placeCells = InfoScoreAll.placeCellID(InfoScoreAll.mouse == i);
    PvaluesAll = [PvaluesAll; Pvalues(placeCells,:)];
    QstatisticsAll = [QstatisticsAll; Qstatistics(placeCells,:)];
end
InfoScoreAll.Pvalues_Ctrl_CNO = PvaluesAll(:,1);
InfoScoreAll.Pvalues_PCtrl_CNO = PvaluesAll(:,2);
InfoScoreAll.Qstatistics_Ctrl_CNO = QstatisticsAll(:,1);
InfoScoreAll.Qstatistics_PCtrl_CNO = QstatisticsAll(:,2);

save placeCellsInfoScoreAllmice_openField_saline_jackknife_nboot5BiningTime.mat InfoScoreAll


%% calculate the percent of each category by varying the threshold of pvalues from jackknife sampling analysis
% load the category information
thresh = 0.05:0.1:1;
percentAll = zeros(length(thresh),3,length(mouse));num_placeCellsAll = zeros(length(thresh),length(mouse));
for ii = 1:length(mouse)
    for i = 1:length(thresh)
        placeCellsSig = find(InfoScoreAll{InfoScoreAll.mouse == ii,15} < thresh(i) & InfoScoreAll{InfoScoreAll.mouse == ii,16} < thresh(i));
        group_new_ii = group_new(InfoScoreAll.mouse == ii);
        num_placeCellsAll(i,ii) = length(placeCellsSig)/length(group_new_ii);
        
        t = group_new_ii(placeCellsSig); t_uni = unique(t);
        for jj = 1:length(t_uni)
            %             percentAll(i,t_uni(jj),ii) = sum(t == t_uni(jj))/length(t);
            percentAll(i,t_uni(jj),ii) = sum(t == t_uni(jj))/length(group_new_ii);
        end
    end
end
percent = mean(percentAll,3);num_placeCells = mean(num_placeCellsAll,2);

colors = [160 255 160;255 160 160;212 212 212;0 0 0]/255;
data = [percent, num_placeCells];
figure
clf
yyaxis left
h = plot(percent,'-o','lineWidth',1.5);
set(h, {'color'}, {colors(1,:); colors(2,:); colors(3,:)});
% ylabel({'% to total number of place cells','passing the significance test'}','FontSize',10,'FontName','Arial')
ylabel('% to total number of place cells','FontSize',10,'FontName','Arial')
[~,h2] = legend(h,{'Bit Decrease','Bit Increase','Un-recovered',''},'Interpreter','none','Location','westoutside','Box','off','FontSize', 10);
set(findobj(h2,'type','patch'),'MarkerSize',10);

yyaxis right
h = plot(num_placeCells,'--o','lineWidth',1.5);
ylabel({'% place cells passing the significance test'}','FontSize',10,'FontName','Arial')
set(gca,'FontSize',10)
xlim([0 length(thresh)+0.5])
xlabels = cellstr(num2str(thresh'));
set(gca,'Xtick',1:1:length(thresh))
set(gca,'XtickLabel',xlabels(1:1:end),'FontSize',10,'FontName','Arial')
xtickangle(30)

xlabel('Threshold of Pvalues (CNO vs Ctrl & Pctrl)', 'FontSize',10,'FontName','Arial')
title('Open field-saline (jackknife, 5 time bins)')

%% bootstrap analysis
load('/Users/suoqinjin/Documents/NN_revision2/Yanjun_openField_results/new_individualShuffle_yanjun0515/placeCellsInfoScoreAllmice_CNOsec_final.mat')
load('/Users/suoqinjin/Documents/NN_revision2/Yanjun_openField_results/new_individualShuffle_yanjun0515/placeCellsInfoScoreAllmice_salinesec_final.mat')
P12All = [];P32All = [];peakFRAll = [];
for i = 1:length(mouse)
    filefolder = fullfile(pwd,mouse{i})
    load(fullfile(filefolder,[mouse{i},'_distributionTest_bootstrap_nboot1000_10BiningTime.mat']),'P12','P32')
    placeCells = InfoScoreAll.placeCellID(InfoScoreAll.mouse == i);
    P12All = [P12All; P12(placeCells,:)];
    P32All = [P32All; P32(placeCells,:)];
end

InfoScoreAll.sig_Ctrl_CNO = P12All;
InfoScoreAll.sig_PCtrl_CNO = P32All;

save placeCellsInfoScoreAllmice_openField_saline_bootstrap_nboot1000_10BiningTime.mat InfoScoreAll


% calculate the percent of each category by varying Confidence Interval of bootstrap analysis
% load the category information
thresh = 0.05:0.1:1;
percentAll = zeros(length(thresh),3,length(mouse));num_placeCellsAll = zeros(length(thresh),length(mouse));
for ii = 1:length(mouse)
    filefolder = fullfile(pwd,mouse{ii})
    load(fullfile(filefolder,[mouse{ii},'_distributionTest_bootstrap_nboot1000_10BiningTime.mat']),'D12','D32')
    for i = 1:length(thresh)
        Q12 = quantile(D12,[thresh(i)/2 1-thresh(i)/2],2);P12 = Q12(:,1).*Q12(:,2) > 0;
        Q32 = quantile(D32,[thresh(i)/2 1-thresh(i)/2],2);P32 = Q32(:,1).*Q32(:,2) > 0;
        placeCells = InfoScoreAll.placeCellID(InfoScoreAll.mouse == ii);
        P12 = P12(placeCells,:);
        P32 = P32(placeCells,:);
        
        placeCellsSig = find(P12 & P32);
        num_placeCellsAll(i,ii) = length(placeCellsSig)/length(placeCells);
        group_new_ii = group_new(InfoScoreAll.mouse == ii);
        t = group_new_ii(placeCellsSig); t_uni = unique(t);
        for jj = 1:length(t_uni)
            %             percentAll(i,t_uni(jj),ii) = sum(t == t_uni(jj))/length(t);
            percentAll(i,t_uni(jj),ii) = sum(t == t_uni(jj))/length(group_new_ii);
        end
    end
end
percent = mean(percentAll,3);num_placeCells = mean(num_placeCellsAll,2);
colors = [160 255 160;255 160 160;212 212 212]/255;
figure
clf
yyaxis left
h = plot(percent,'-o','lineWidth',1.5);
set(h, {'color'}, {colors(1,:); colors(2,:); colors(3,:)});
% ylabel({'% to total number of place cells','passing the significance test'}','FontSize',10,'FontName','Arial')
ylabel('% to total number of place cells','FontSize',10,'FontName','Arial')
[~,h2] = legend(h,{'Bit Decrease','Bit Increase','Un-recovered',''},'Interpreter','none','Location','westoutside','Box','off','FontSize', 10);
set(findobj(h2,'type','patch'),'MarkerSize',10);

yyaxis right
h = plot(num_placeCells,'--o','lineWidth',1.5);
ylabel({'% place cells passing the significance test'}','FontSize',10,'FontName','Arial')
set(gca,'FontSize',10)
xlim([0 length(thresh)+0.5])
xlabels = cellstr(num2str(thresh'));
set(gca,'Xtick',1:1:length(thresh))
set(gca,'XtickLabel',xlabels(1:1:end),'FontSize',10,'FontName','Arial')
xtickangle(30)

xlabel('100*(1-alpha) bootstrap confidence interval (CNO vs Ctrl & Pctrl)', 'FontSize',10,'FontName','Arial')
title('Open field-saline (1000 bootstrap, 10 time bins)')



