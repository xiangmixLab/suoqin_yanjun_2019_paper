%%%%%%%%% Calcium imaging data analysis (April-2019)-(updated July-2019) %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% README start %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code load InfoScoreAll variable from place field analysis, please see the general pipelin ein place_cell_analysis_OLM_basedTrainingPCs.m
% Functionality of this code: 
% (1) compare the spatial information (information score) training-testing in saline vs training-testing in CNO (population level)
% (2) classifying place cells based on the response to object
% (3) comparison of the difference in percent between obj1 and obj2, between saline and CNO 
% (4) compare the in-field peak rate
%%%%%%%%%%%%%%%%%%%% README end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load the data
mice_saline = {'M3421F','M3423F','M3427F','M3411','M3414'};
mice_CNO = {'M3425F','M3426F','M3422F','M3424F','M3415','M3412'};
mice_rev_saline = {'M3425F_rev','M3413_rev','M3424F_rev','M3415_rev','M3412_rev'};
mice_rev_CNO = {'M3421F_rev','M3423F_rev','M3422F_rev','M3411_rev','M3414_rev'};
mouse = [mice_saline,mice_rev_saline];
mouse = [mice_CNO,mice_rev_CNO];
sessions = {'Ctrl','CNO','PCtrl'};
%% category analysis based on bootstrapping
load('placeCellsInfoScoreAllmice_OLM_placeFieldAnalysis_individualPlaceCells_testing_saline.mat')
load('placeCellsInfoScoreAllmice_OLM_placeFieldAnalysis_individualPlaceCells_testing_CNO.mat')
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

%% compare the spatial information (information score) between the (Ctrl - Post-Crtl) and (Ctrl - CNO: difference between Ctrl and CNO) and .
data_InfoScoreAll_diff = [data_InfoScoreAll(:,1) - data_InfoScoreAll(:,3),data_InfoScoreAll(:,1) - data_InfoScoreAll(:,2)];
colors = [39 170 225;213 128 43]/255;

data = (data_InfoScoreAll_diff);

figure
clf
h = boxplot(data,'colors',colors,'symbol', '.','OutlierSize',10);
set(gca,'FontSize',8)
xlim([0.5 size(data,2)+0.5])
xlabels = {'Ctrl-PCtrl','Ctrl-CNO'};
set(gca,'Xtick',1:size(data,2))
set(gca,'XtickLabel',xlabels,'FontSize',10,'FontName','Arial')
%xtickangle(30)
ylabel('Bit / Sec','FontSize',10,'FontName','Arial')
% Alter linestyle
idxColor = fliplr(1:size(data,2));
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
% ranksum(data(:,1),data(:,2),'tail','left')
% [~,p] = ttest2(data(:,1),data(:,2),'tail','left')
[~,p] = ttest(data(:,1),data(:,2))
[~,p] = ttest(data(:,1))
[~,p] = ttest(data(:,2))
text(0+0.35,max(data(:))+0.02,['P=' num2str(p,'%3.0e')],'FontSize',10)
text(0+0.35,max(data(:))+0.02,['P=' num2str(p,'%.3f')],'FontSize',10)


%% classifying place cells
colors_cluster = [228,26,28; 55,126,184; 212 212 212]/255;
colors_cluster = [228,26,28; 55,126,184; 152,78,163; 212 212 212]/255;

group = zeros(size(data_InfoScoreAll,1),1);
group(sum(data_InfoScoreAll,2) == 2) = 3;
group(data_InfoScoreAll(:,3) == 1) = 4; 
group(data_InfoScoreAll(:,1) == 1 & group ~= 3) = 1;
group(data_InfoScoreAll(:,2) == 1 & group ~= 3) = 2;

%% calculate the percent of each group
% % percent = zeros(3,length(mouse));
% % for i = 1:length(mouse)
% %      percent(:,i) = mean(data_InfoScoreAll(InfoScoreAll.mouse == i,:));
% % end
% % 
% % percent = zeros(4,length(mouse));
% % for i = 1:length(mouse)
% %     percent(3,i) = sum(sum(data_InfoScoreAll(InfoScoreAll.mouse == i,:),2) == 2)/sum(InfoScoreAll.mouse == i);
% %     percent(1:2,i) = sum(data_InfoScoreAll(InfoScoreAll.mouse == i,1:2))/sum(InfoScoreAll.mouse == i)-percent(3,i);
% %     percent(4,i) = sum(data_InfoScoreAll(InfoScoreAll.mouse == i,3))/sum(InfoScoreAll.mouse == i);
% % end

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

% compare the in-field peak rate
colors = colors_cluster; 
data_PR= [InfoScoreAll{:,8},InfoScoreAll{:,9},InfoScoreAll{:,10}]; %
data = sum(data_PR,2); data(group == 3) = mean(data_PR(group == 3,1:2),2);

data = InfoScoreAll{:,4}; 
figure
clf
h = boxplot(data,group,'colors',colors,'symbol', '.','OutlierSize',10);
xlim([0.5 size(colors,1)+0.5])
xlabels = {'obj 1','obj 2','both','out-field'};
set(gca,'Xtick',1:size(colors,1))
set(gca,'XtickLabel',xlabels,'FontSize',10,'FontName','Arial')
xtickangle(30)
ylabel('Peak in-field event rate','FontSize',10,'FontName','Arial')
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
text(0+0.35,max(data(:))+0.02,['P=' num2str(p,'%3.0e')],'FontSize',10)
text(1+0.35,max(data(:))+0.02,['P=' num2str(p,'%.3f')],'FontSize',10)



% compare the information score between the (Ctrl - Post-Crtl) and (Ctrl - CNO) in different groups
colors = [39 170 225;213 128 43]/255;
colorsT = [colors;0 0 0;colors;0 0 0;colors];
figure
clf
% boxplot([data_InfoScoreAll_diff(InfoScoreAll.group == 1,:),nan(nnz(InfoScoreAll.group == 1),6)],'colors',colorsT,'symbol', '.')
% hold on
% boxplot([nan(nnz(InfoScoreAll.group == 2),3),data_InfoScoreAll_diff(InfoScoreAll.group == 2,:),nan(nnz(InfoScoreAll.group == 2),3)],'colors',colorsT,'symbol', '.')
% boxplot([nan(nnz(InfoScoreAll.group == 3),6),data_InfoScoreAll_diff(InfoScoreAll.group == 3,:)],'colors',colorsT,'symbol', '.')
boxplot([data_InfoScoreAll_diff(group == 1,:),nan(nnz(group == 1),6)],'colors',colorsT,'symbol', '.','OutlierSize',10)
hold on
boxplot([nan(nnz(group == 2),3),data_InfoScoreAll_diff(group == 2,:),nan(nnz(group == 2),3)],'colors',colorsT,'symbol', '.','OutlierSize',10)
boxplot([nan(nnz(group == 3),6),data_InfoScoreAll_diff(group == 3,:)],'colors',colorsT,'symbol', '.','OutlierSize',10)
%legend(h([1 2]),{'Ctrl-PCtrl','Ctrl-CNO'},'FontSize',10)
set(gca,'FontSize',8)
xlim([0.5 8+0.5])
xlabels = {'Bit Decrease','Bit Increase','Un-recovered'};
set(gca,'Xtick',[1.5 4.5 7.5])
set(gca,'XtickLabel',xlabels,'FontSize',10,'FontName','Arial')
% xtickangle(30)
ylabel('Bit / Sec','FontSize',10,'FontName','Arial')
ylabel('Ca^{2+} event rate (Hz)','FontSize',10,'FontName','Arial')
% Alter linestyle
colorsTT = [colorsT;colorsT;colorsT];
idxColor = fliplr(1:size(colorsTT,1));
h2 = findobj(gca,'Tag','Box');
for j=1:length(h2)
    patch(get(h2(j),'XData'),get(h2(j),'YData'),colorsTT(idxColor(j),:),'FaceAlpha',0.6);
end
h3 = findobj(gca,'tag','Outliers');
for j = 1:length(h3)
    h3(j).MarkerEdgeColor = colorsTT(idxColor(j),:);
end
h1 = findobj(gca,'tag','Median');
set(h1,{'linew'},{2.5})
set(h1,'Color',[0 0 0])
box off


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

% xtickangle(30)
[~,hle] = legend(h2([2 1]),{'Ctrl-PCtrl','Ctrl-CNO'},'FontSize',10);
hl = findobj(hle,'type','line');
set(hl,'LineWidth',8);

%[~,p] = ttest(data_InfoScoreAll_diff(group_new == 1,1),data_InfoScoreAll_diff(group_new == 1,2),'tail','left')
[~,p] = ttest(data_InfoScoreAll_diff(group == 1,1),data_InfoScoreAll_diff(group == 1,2)) % change to two-tailed in final version 03-03-2019
%text(1+0.35,1+0.02,['P=' num2str(p,'%.2f')],'FontSize',10)
text(1+0.35,1+0.02,['P=' num2str(p,'%3.0e')],'FontSize',8)
%[~,p] = ttest(data_InfoScoreAll_diff(group_new == 2,1),data_InfoScoreAll_diff(group_new == 2,2),'tail','right')
[~,p] = ttest(data_InfoScoreAll_diff(group == 2,1),data_InfoScoreAll_diff(group == 2,2))
%text(4+0.35,1+0.02,['P=' num2str(p,'%.2f')],'FontSize',10)
text(4+0.35,1+0.02,['P=' num2str(p,'%3.0e')],'FontSize',8)
%[~,p] = ttest(data_InfoScoreAll_diff(group_new == 3,1),data_InfoScoreAll_diff(group_new == 3,2),'tail','right')
[~,p] = ttest(data_InfoScoreAll_diff(group == 3,1),data_InfoScoreAll_diff(group == 3,2))
text(7+0.35,1+0.02,['P=' num2str(p,'%.2f')],'FontSize',8)

%% compare the information score between the (Ctrl - Post-Crtl) and (Ctrl - CNO) in different groups
% data_InfoScoreAll = [InfoScoreAll{:,5},InfoScoreAll{:,6},InfoScoreAll{:,7}]; % FR
data_InfoScoreAll = [InfoScoreAll{:,2},InfoScoreAll{:,3},InfoScoreAll{:,4}]; % SI
data_InfoScoreAll_diff = [data_InfoScoreAll(:,1) - data_InfoScoreAll(:,3),data_InfoScoreAll(:,1) - data_InfoScoreAll(:,2)];


