%%%%%%%%% Calcium imaging data analysis (April-2019)-(updated July-2019) %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% README start %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code is a general pipeline for place field analysis in OLM experiment, which load InfoScoreAll variable from place field analysis
% Functionality of this code: 
% (1) compare the spatial information (information score) training-testing in saline vs training-testing in CNO (population level)
% (2) classifying place cells based on the response to object
% (3) comparison of the difference in percent between obj1 and obj2, between saline and CNO 
% (4)  compare the in-field peak rate
%%%%%%%%%%%%%%%%%%%% README end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load the data
mice_saline = {'M3421F','M3423F','M3427F','M3411','M3414'};
mice_CNO = {'M3425F','M3426F','M3422F','M3424F','M3415','M3412'};
mice_rev_saline = {'M3425F_rev','M3413_rev','M3424F_rev','M3415_rev','M3412_rev'};
mice_rev_CNO = {'M3421F_rev','M3423F_rev','M3422F_rev','M3411_rev','M3414_rev'};
mouse = [mice_saline,mice_rev_saline];
mouse = [mice_CNO,mice_rev_CNO];
sessions = {'Ctrl','CNO','PCtrl'};

load('placeCellsInfoScoreAllmice_OLM_placeFieldAnalysis_basedTrainingPC_saline.mat')
load('placeCellsInfoScoreAllmice_OLM_placeFieldAnalysis_basedTrainingPC_CNO.mat')
%% get the calculated info scores
data_PR1 = [InfoScoreAll{:,14},InfoScoreAll{:,15},InfoScoreAll{:,16}]; % (training)In-field peak event rate, obj1/obj2/out-field
data_PR2 = [InfoScoreAll{:,17},InfoScoreAll{:,18},InfoScoreAll{:,19}]; % (testing)In-field peak event rate, obj1/obj2/out-field

% filter neurons
InfoScoreAll0 = InfoScoreAll; 
neuron_lowRate1 = find(sum(data_PR1 < 0.1,2) == 3);
neuron_lowRate2 = find(sum(data_PR2 < 0.1,2) == 3);
neuron_lowRate = intersect(neuron_lowRate1, neuron_lowRate2);

%neuron_lowRate = find(InfoScoreAll{:,4} < 0.01);
InfoScoreAll(neuron_lowRate,:) = [];

data_InfoScoreAll = [InfoScoreAll{:,2},InfoScoreAll{:,3},InfoScoreAll{:,4}]; % SI-bits/sec
data_InfoScoreAll = [InfoScoreAll{:,5},InfoScoreAll{:,6},InfoScoreAll{:,7}]; % SI-bits/spike
data_InfoScoreAll = [InfoScoreAll{:,8},InfoScoreAll{:,9},InfoScoreAll{:,10}]; % training place field response, obj1/obj2/out-field
data_InfoScoreAll = [InfoScoreAll{:,11},InfoScoreAll{:,12},InfoScoreAll{:,13}]; % testing place field response, obj1/obj2/out-field

%% compare the spatial information (information score) training-testing in saline vs training-testing in CNO (population level)

InfoScoreAllAll_CNO = [ones(size(InfoScoreAll_CNO,1),1), InfoScoreAll_CNO{:,2:3}]; % bit/sec
InfoScoreAllAll_Saline = [ones(size(InfoScoreAll_saline,1),1), InfoScoreAll_saline{:,2:3}]; % bits/sec
InfoScoreAllAll_CNO = [ones(size(InfoScoreAll_CNO,1),1), InfoScoreAll_CNO{:,4:5}]; % bits/spike
InfoScoreAllAll_Saline = [ones(size(InfoScoreAll_saline,1),1), InfoScoreAll_saline{:,4:5}]; % bits/spike
MeanFiringRateAllAll_CNO = [ones(size(InfoScoreAll_CNO,1),1), InfoScoreAll_CNO{:,6:7}]; % mean ER
MeanFiringRateAllAll_Saline = [ones(size(InfoScoreAll_saline,1),1), InfoScoreAll_saline{:,6:7}]; % mean ER
conditions = {'Saline','CNO'};show_pvalue = 1; measure = 'SI';
colors_conditions = [228,26,28;55,126,184]/255;
SpatialInfoComparison_conditions(InfoScoreAllAll_Saline,InfoScoreAllAll_CNO,colors_conditions,conditions,show_pvalue,measure)
SpatialInfoComparison_conditions(MeanFiringRateAllAll_Saline,MeanFiringRateAllAll_CNO,colors_conditions,conditions,show_pvalue,'event rate')

data_InfoScoreAll_diff = [data_InfoScoreAll(:,1) - data_InfoScoreAll(:,3),data_InfoScoreAll(:,1) - data_InfoScoreAll(:,2)];


%% classifying place cells based on the response to object
data_InfoScoreAll = [InfoScoreAll{:,8},InfoScoreAll{:,9},InfoScoreAll{:,10}]; % training place field response, obj1/obj2/out-field
data_InfoScoreAll = [InfoScoreAll{:,11},InfoScoreAll{:,12},InfoScoreAll{:,13}]; % testing place field response, obj1/obj2/out-field
colors_cluster = [228,26,28; 55,126,184; 212 212 212]/255;
colors_cluster = [228,26,28; 55,126,184; 152,78,163; 212 212 212]/255;

% 1) response to object 1; 2) response to object 2; 3) response to both objects; 4) response to other places
group = zeros(size(data_InfoScoreAll,1),1);
group(sum(data_InfoScoreAll,2) == 2) = 3;
group(data_InfoScoreAll(:,3) == 1) = 4; 
group(data_InfoScoreAll(:,1) == 1 & group ~= 3) = 1;
group(data_InfoScoreAll(:,2) == 1 & group ~= 3) = 2;

% 1) response to object 1; 2) response to object 2; 3) response to other places
group = zeros(size(data_InfoScoreAll,1),1);
group(data_InfoScoreAll(:,3) == 1) = 3; 
group(data_InfoScoreAll(:,1) == 1) = 1;
group(data_InfoScoreAll(:,2) == 1) = 2;

%% calculate the percent of each category group
percent = zeros(4,length(mouse));
for i = 1:length(mouse)
    group_ii = group(InfoScoreAll.mouse == i);
    t = group_ii; t_uni = unique(t);
    for jj = 1:length(t_uni)
        percent(t_uni(jj),i) = sum(t == t_uni(jj))/length(group_ii);
    end
end

mean(percent,2) % 0.225, 0.258, 0.59 saline testing
mean(percent,2) % 0.152, 0.185, 0.07,0.59 saline testing
mean(percent,2) % 0.18, 0.15, 0.06, 0.61 CNO testing

InfoScoreAll_CNO = InfoScoreAll;
group_testing_CNO = group;
percent_testing_CNO = percent;


data = percent'*100;


% bar plot
figure
clf
hBar = bar(mean(data));
ctr = [];ydt = [];
for k1 = 1:length(hBar)
    ctr(k1,:) = bsxfun(@plus, hBar(1).XData, [hBar(k1).XOffset]');
    ydt(k1,:) = hBar(k1).YData;
end
hBar.FaceColor = 'flat';
for j = 1:size(hBar.CData,1)
    hBar.CData(j,:) = colors_cluster(j,:);
end
xlim([0.5 size(hBar.CData,1)+0.5])
plotSpread(data,'distributionMarkers', {'o', 'o', 'o','o'},'distributionColors','k')

xlim([0.5 size(data,2)+0.5])
xlabels = {'obj 1','obj 2','out-field'};
xlabels = {'obj 1','obj 2','both','out-field'};
set(gca,'Xtick',1:size(data,2))
set(gca,'XtickLabel',xlabels,'FontSize',12,'FontName','Arial')
xtickangle(30)
ylabel({'% of total place cells'},'FontSize',12,'FontName','Arial')
ytickformat('%g%%');
set(gca,'FontSize',12,'FontName','Arial')
hold on
errorbar(ctr, ydt, std(data), '.k','marker', 'none','LineWidth',1)
hold off
box off
set(gca,'linewidth',1.5)
[~,p] = ttest(percent(1,:),percent(2,:))
text(1,60,['P=' num2str(p,'%.2f')],'FontSize',8)
%text(1,0.6,['P=' num2str(p,'%3.0e')],'FontSize',8)
title('Saline (Testing)', 'FontSize',12,'FontName','Arial')


%% comparison of the difference in percent between obj1 and obj2, between saline and CNO 
data = [percent_testing(1,:)-percent_testing(2,:), percent_testing_CNO(1,:)-percent_testing_CNO(2,:)];
data = [percent_training(1,:)-percent_training(2,:), percent_training_CNO(1,:)-percent_training_CNO(2,:)];
group = ones(length(data),1)*2; group(1:length(percent_testing(1,:))) = 1;

colors_conditions = [228,26,28;55,126,184]/255;conditions = {'Saline','CNO'};show_pvalue = 1; 
hFig = figure('position', [200, 200, 140,160]);
h = boxplot(data,group,'symbol', '.','OutlierSize',10);

xlim([0.5 2+0.5])
set(gca,'Xtick',1:2)
set(gca,'XtickLabel',conditions,'FontSize',12,'FontName','Arial')
xtickangle(0)

ylabel({'Difference in the percent of place', 'cells between obj 1 and obj 2'},'FontSize',12,'FontName','Arial')
set(gca,'FontSize',12,'FontName','Arial')
% Alter linestyle
idxColor = fliplr(1:2);
h2 = findobj(gca,'Tag','Box');
for jj=1:length(h2)
    patch(get(h2(jj),'XData'),get(h2(jj),'YData'),colors_conditions(idxColor(jj),:),'FaceAlpha',0.6);
end
h3 = findobj(gca,'tag','Outliers');
for jj = 1:length(h3)
    h3(jj).MarkerEdgeColor = colors_conditions(idxColor(jj),:);
end
box off
lines = findobj(gcf,'Type','Line');
set(lines,'LineWidth',1)
h1 = findobj(gcf,'tag','Median');
set(h1,{'linew'},{2.5})
set(h1,'Color',[0 0 0])
h = findall(gcf,'marker','.');
set(h, 'markersize', 10)

set(gca,'LineWidth',1)
set(gca,'FontName','Arial')
set(gca,'FontSize',12)

txt_obj = findall(gca,'Type','text');
set(txt_obj,'FontName','Arial','FontSize',12);

hold on
if show_pvalue
    for i = 1:2-1
        p = ranksum(data(group == 1),data(group == 2))
        p = anova1(data,group, 'off')
        if p < 0.01
            text(i+0.1,max(data(:))+0.01,['P=' num2str(p,'%3.0e')],'FontSize',8,'FontWeight','bold')
        else
            text(i+0.1,max(data(:))+0.01,['P=' num2str(p,'%.3f')],'FontSize',8,'FontWeight','bold')
        end
    end
end


%% compare the in-field peak rate
colors = colors_cluster; 
data_PR1 = [InfoScoreAll{:,14},InfoScoreAll{:,15},InfoScoreAll{:,16}]; % (training)In-field peak event rate, obj1/obj2/out-field
data_PR2 = [InfoScoreAll{:,17},InfoScoreAll{:,18},InfoScoreAll{:,19}]; % (testing)In-field peak event rate, obj1/obj2/out-field
data_PR = data_PR1; 
data = sum(data_PR,2); data(group == 3) = mean(data_PR(group == 3,1:2),2);

data = InfoScoreAll_saline{:,21}; % peak event rate

figure
clf
h = boxplot(data,group,'colors',colors,'symbol', '.','OutlierSize',10);
xlim([0.5 size(colors,1)+0.5])
xlabels = {'obj 1','obj 2','both','out-field'};
set(gca,'Xtick',1:size(colors,1))
set(gca,'XtickLabel',xlabels,'FontSize',10,'FontName','Arial')
xtickangle(30)
ylabel('Peak event rate','FontSize',10,'FontName','Arial')
title('CNO (Testing)','FontSize',10,'FontName','Arial')
% Alter linestyle
idxColor = fliplr(1:size(colors,1));
h2 = findobj(gca,'Tag','Box');
for j=1:length(h2)
    patch(get(h2(j),'XData'),get(h2(j),'YData'),colors(idxColor(j),:),'FaceAlpha',0.6);
end
h3 = findobj(gca,'tag','Outliers');
for j = 1:length(h3)
    h3(j).MarkerEdgeColor = colors(idxColor(j),:);
end
set(h,'LineWidth',1)
set(gca,'FontSize',12)
set(gca,'linewidth',1.5)

h1 = findobj(gca,'tag','Median');
set(h1,{'linew'},{2.5})
set(h1,'Color',[0 0 0])
box off

h1 = findobj(gca,'tag','Upper Whisker');
set(h1,'LineWidth',1)
set(h1,'LineStyle','-')
h1 = findobj(gca,'tag','Lower Whisker');
set(h1,'LineWidth',1)
set(h1,'LineStyle','-')
p = ranksum(data(group == 1),data(group == 2))
%text(0+0.35,max(data(:))+0.02,['P=' num2str(p,'%3.0e')],'FontSize',10)
text(1+0.35,max(data(:))+0.02,['P=' num2str(p,'%.3f')],'FontSize',10)


