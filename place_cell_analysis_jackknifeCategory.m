%%%%%%%%% Calcium imaging data analysis (Jan-2019)-(April-2019)-(updated July-2019) %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% README start %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code is a general pipeline for category analysis of place cells, which load InfoScoreAll variable
% Functionality of this code:
% (1) compare the spatial information (information score) between the (Ctrl - Post-Crtl) and (Ctrl - CNO: difference between Ctrl and CNO)
% (2) classifying place cells based on the change of SI and jackknife analysis
% (3) compare the information score or event rate between the (Ctrl - Post-Crtl) and (Ctrl - CNO) in different groups
%%%%%%%%%%%%%%%%%%%% README end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% load the data
mouse = {'M3321','1stCA1','M3322'}; % saline mice
mouse = {'M3244F','M3321','1stCA1','M3322','M3243','M3323b'}; % CNO mice
sessions = {'Ctrl','CNO','PCtrl'};
%% category analysis based on bootstrapping
% load('placeCellsInfoScoreAllmice_linearTrack_CNO_bootstrap_nboot1000_10BiningTime_bitSpikeYanjun.mat')
load('placeCellsInfoScoreAllmice_linearTrack_CNO_jackknife_nboot10BiningTime_bitSecYanjun.mat')
load('placeCellsInfoScoreAllmice_linearTrack_saline_jackknife_nboot10BiningTime_bitSecYanjun.mat')
%% get the calculated info scores
data_InfoScoreAll = [InfoScoreAll{:,2},InfoScoreAll{:,3},InfoScoreAll{:,4}]; % SI-bits/sec
data_InfoScoreAll = [InfoScoreAll{:,5},InfoScoreAll{:,6},InfoScoreAll{:,7}]; % SI-bits/spike
data_InfoScoreAll = [InfoScoreAll{:,11},InfoScoreAll{:,12},InfoScoreAll{:,13}]; % Peak event rate
data_InfoScoreAll = [InfoScoreAll{:,8},InfoScoreAll{:,9},InfoScoreAll{:,10}]; % Mean event rate
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


%% classifying place cells based on the change of SI and jackknife analysis
colors_cluster = [160 255 160;255 160 160;166,206,227;212 212 212]/255;
% define a new group based on the scatter plot and jackknife based t-test
group_new = ones(size(data_InfoScoreAll_diff,1),1)*3;
group_new(data_InfoScoreAll_diff(:,2) > 0 & data_InfoScoreAll_diff(:,2) > data_InfoScoreAll_diff(:,1)) = 1;
group_new(data_InfoScoreAll_diff(:,2) < 0 & data_InfoScoreAll_diff(:,2) < data_InfoScoreAll_diff(:,1)) = 2;

%  group_new(InfoScoreAll{:,15} == 0 | InfoScoreAll{:,16} ==0) = 4;
group_new(InfoScoreAll{:,15} > 0.05 | InfoScoreAll{:,16} > 0.05) = 4; % un-assigned group

figure
h = [];
for i = 1:size(colors_cluster,1)
    h(i) = scatter(data_InfoScoreAll_diff(group_new == i,1),data_InfoScoreAll_diff(group_new == i,2),10,colors_cluster(i,:),'filled');
    hold on
    line([-2 2],[-2 2],'LineStyle','--','color','k','LineWidth',1)
    line([-2 2],[0 0],'LineStyle','--','color','k','LineWidth',1)
end
set(gca,'FontSize',12)
xlabel('Ctrl-PCtrl (Bit / Sec)','FontSize',12,'FontName','Arial');
ylabel('Ctrl-CNO (Bit / Sec)','FontSize',12,'FontName','Arial');
%axis([min(ydata(:,1))-0.01 max(ydata(:,1))+0.01, min(ydata(:,2))-0.01 max(ydata(:,2))+0.01])
set(gca,'linewidth',1.5)
box off
[~,h2] = legend(h,{'Bit Decrease','Bit Increase','Un-recovered','Un-assigned'},'Interpreter','none','Location','eastoutside','Box','off','FontSize', 10);
set(findobj(h2,'type','patch'),'MarkerSize',10);

%% calculate percent in each group
x = data_InfoScoreAll_diff(:,1); y = data_InfoScoreAll_diff(:,2);
q1 = x > 0 & y > 0 & x < y ; q2 = x > 0 & y > 0 & x > y; q3 = x > 0 & y < 0 & x > -y; q4 = x > 0 & y < 0 & x < -y;
q5 = x < 0 & x > y & y < 0; q6 = x < 0 & x < y & y < 0; q7 = y > 0 & x < -y & x < 0; q8 = y > 0 & x > -y & x < 0;
group_p = zeros(size(data_InfoScoreAll_diff,1),1);
group_p(q1) = 1; group_p(q2) = 2; group_p(q3) = 3; group_p(q4) = 4;
group_p(q5) = 5; group_p(q6) = 6; group_p(q7) = 7; group_p(q8) = 8;
percent_p = grpstats(group_p,group_p,'numel')/length(group_p);
text_pos = [0.5 1.5; 1.5 0.5; 1.5 -0.5; 0.5 -1.5; -0.5 -1.5; -1.5 -0.5; -1.5 0.5; -0.5 1.5];
for i = 1:length(percent_p)
    text(text_pos(i,1),text_pos(i,2),num2str(round(percent_p (i)*100),'%g%%'),'FontSize',10,'FontName','Arial')
    hold on
end

figure
h = [];
for i = 1:size(colors_cluster,1)
    h(i) = scatter(data_InfoScoreAll_diff(group_new == i,1),data_InfoScoreAll_diff(group_new == i,2),15,colors_cluster(i,:),'filled');
    hold on
    line([-2 2],[-2 2],'LineStyle','--','color','k')
    line([-2 2],[2 -2],'LineStyle','--','color','k')
    line([-2 2],[-2 2],'LineStyle','--','color','k')
    line([-2 2],[0 0],'LineStyle','--','color','k')
    line([0 0],[-2 2],'LineStyle','--','color','k')
end
for i = 1:length(percent_p)
    text(text_pos(i,1),text_pos(i,2),num2str(round(percent_p (i)*100),'%g%%'),'FontSize',10,'FontName','Arial')
    hold on
end
set(gca,'FontSize',8)
% xlabel('Ctrl-PCtrl (Bit / Sec)','FontSize',10,'FontName','Arial');
% ylabel('Ctrl-CNO (Bit / Sec)','FontSize',10,'FontName','Arial');
xlabel('CNO-Ctrl (Bit / Sec)','FontSize',10,'FontName','Arial');
ylabel('CNO-PCtrl (Bit / Sec)','FontSize',10,'FontName','Arial');
title('Neurons are colored by new category','FontSize',10)
box off
[~,h2] = legend(h,{'Bit Decrease','Bit Increase','Un-recovered'},'Interpreter','none','Location','eastoutside','Box','off','FontSize', 10);
set(findobj(h2,'type','patch'),'MarkerSize',10);



%% calculate the percent of decreased/increased group
% percent = zeros(3,length(mouse));
% for i = 1:length(mouse)
%     %     persent(:,i) = grpstats(InfoScore{i}.group,InfoScore{i}.group,'numel')/length(InfoScore{i}.group);
%     percent(:,i) = grpstats(group_new(InfoScoreAll.mouse == i),group_new(InfoScoreAll.mouse == i),'numel')/length(group_new(InfoScoreAll.mouse == i));
% end
percent = zeros(3,length(mouse));
for i = 1:length(mouse)
    % persent(:,i) = grpstats(InfoScore{i}.group,InfoScore{i}.group,'numel')/length(InfoScore{i}.group);
    % percent(:,i) = grpstats(group_new(InfoScoreAll.mouse == i),group_new(InfoScoreAll.mouse == i),'numel')/length(group_new(InfoScoreAll.mouse == i));
    group_new_ii = group_new(InfoScoreAll.mouse == i);
    t = group_new_ii; t_uni = unique(t);
    for jj = 1:length(t_uni)
        percent(t_uni(jj),i) = sum(t == t_uni(jj))/length(group_new_ii);
    end
end
mean(percent,2) % 0.4566 0.2565 0.2869
mean(percent,2) % 0.4135 0.2433 0.3432
% percent = grpstats(group_new,group_new,'numel')/length(group_new);
data = percent'*100;
% box plot
figure
clf
h = boxplot(data,'colors',colors,'symbol', '.');
set(gca,'FontSize',8)
xlim([0.5 size(data,2)+0.5])
xlabels = {'Bit Decrease','Bit Increase','Un-recovered'};
set(gca,'Xtick',1:size(data,2))
set(gca,'XtickLabel',xlabels,'FontSize',10,'FontName','Arial')
xtickangle(30)
ylabel('% to total number of place cells','FontSize',10,'FontName','Arial')
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
h1 = findobj(gca,'tag','Median');
set(h1,{'linew'},{2.5})
set(h1,'Color',[0 0 0])
box off

[~,p] = ttest(percent(1,:),percent(2,:),'tail','right')
%text(1+0.35,1+0.02,['P=' num2str(p,'%.2f')],'FontSize',10)
text(1+0.35,1+0.02,['P=' num2str(p,'%3.0e')],'FontSize',8)

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
plotSpread(data,'distributionMarkers', {'o', 'o', 'o'},'distributionColors','k')

xlim([0.5 size(data,2)+0.5])
% xlabels = {'Bit Decrease','Bit Increase','Un-recovered'};
xlabels = {'Bit Decrease','Bit Increase','Un-recovered','Un-assigned'};
set(gca,'Xtick',1:size(data,2))
set(gca,'XtickLabel',xlabels,'FontSize',12,'FontName','Arial')
xtickangle(30)
ylabel({'% of place cells showing', 'significant difference'},'FontSize',12,'FontName','Arial')
ytickformat('%g%%');
set(gca,'FontSize',12,'FontName','Arial')
hold on
errorbar(ctr, ydt, std(data), '.k','marker', 'none','LineWidth',1)
hold off
box off
set(gca,'linewidth',1.5)
[~,p] = ttest(percent(1,:),percent(2,:))
text(1,0.6,['P=' num2str(p,'%.2f')],'FontSize',8)
%text(1,0.6,['P=' num2str(p,'%3.0e')],'FontSize',8)

% compare the information score between the (Ctrl - Post-Crtl) and (Ctrl - CNO) in different groups
colors = [39 170 225;213 128 43]/255;
colorsT = [colors;0 0 0;colors;0 0 0;colors];
figure
clf
% boxplot([data_InfoScoreAll_diff(InfoScoreAll.group == 1,:),nan(nnz(InfoScoreAll.group == 1),6)],'colors',colorsT,'symbol', '.')
% hold on
% boxplot([nan(nnz(InfoScoreAll.group == 2),3),data_InfoScoreAll_diff(InfoScoreAll.group == 2,:),nan(nnz(InfoScoreAll.group == 2),3)],'colors',colorsT,'symbol', '.')
% boxplot([nan(nnz(InfoScoreAll.group == 3),6),data_InfoScoreAll_diff(InfoScoreAll.group == 3,:)],'colors',colorsT,'symbol', '.')
boxplot([data_InfoScoreAll_diff(group_new == 1,:),nan(nnz(group_new == 1),6)],'colors',colorsT,'symbol', '.','OutlierSize',10)
hold on
boxplot([nan(nnz(group_new == 2),3),data_InfoScoreAll_diff(group_new == 2,:),nan(nnz(group_new == 2),3)],'colors',colorsT,'symbol', '.','OutlierSize',10)
boxplot([nan(nnz(group_new == 3),6),data_InfoScoreAll_diff(group_new == 3,:)],'colors',colorsT,'symbol', '.','OutlierSize',10)
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
[~,p] = ttest(data_InfoScoreAll_diff(group_new == 1,1),data_InfoScoreAll_diff(group_new == 1,2)) % change to two-tailed in final version 03-03-2019
%text(1+0.35,1+0.02,['P=' num2str(p,'%.2f')],'FontSize',10)
text(1+0.35,1+0.02,['P=' num2str(p,'%3.0e')],'FontSize',8)
%[~,p] = ttest(data_InfoScoreAll_diff(group_new == 2,1),data_InfoScoreAll_diff(group_new == 2,2),'tail','right')
[~,p] = ttest(data_InfoScoreAll_diff(group_new == 2,1),data_InfoScoreAll_diff(group_new == 2,2))
%text(4+0.35,1+0.02,['P=' num2str(p,'%.2f')],'FontSize',10)
text(4+0.35,1+0.02,['P=' num2str(p,'%3.0e')],'FontSize',8)
%[~,p] = ttest(data_InfoScoreAll_diff(group_new == 3,1),data_InfoScoreAll_diff(group_new == 3,2),'tail','right')
[~,p] = ttest(data_InfoScoreAll_diff(group_new == 3,1),data_InfoScoreAll_diff(group_new == 3,2))
text(7+0.35,1+0.02,['P=' num2str(p,'%.2f')],'FontSize',8)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The followings are not used in the paper

%% compare the information score between the (Ctrl - Post-Crtl) and (Ctrl - CNO) in different groups
% data_InfoScoreAll = [InfoScoreAll{:,5},InfoScoreAll{:,6},InfoScoreAll{:,7}]; % FR
data_InfoScoreAll = [InfoScoreAll{:,2},InfoScoreAll{:,3},InfoScoreAll{:,4}]; % SI
data_InfoScoreAll_diff = [data_InfoScoreAll(:,1) - data_InfoScoreAll(:,3),data_InfoScoreAll(:,1) - data_InfoScoreAll(:,2)];

%% compute the raio shift
data_InfoScoreAll_diff = [data_InfoScoreAll(:,1) - data_InfoScoreAll(:,3),data_InfoScoreAll(:,1) - data_InfoScoreAll(:,2)];
data_InfoScoreAll_diff_ratio = data_InfoScoreAll_diff(:,2)./data_InfoScoreAll_diff(:,1)-1;

data_InfoScoreAll_diffPercent = [data_InfoScoreAll(:,3)./data_InfoScoreAll(:,1) - 1, data_InfoScoreAll(:,2)./data_InfoScoreAll(:,1) - 1];
data_InfoScoreAll_diffPercent_relative = data_InfoScoreAll_diffPercent(:,2)-data_InfoScoreAll_diffPercent(:,1);

data = abs(data_InfoScoreAll_diffPercent_relative)*100;
idx = data_InfoScoreAll_diffPercent(:,1).*data_InfoScoreAll_diffPercent(:,2) > 0;
data = data(idx);
figure
boxplot(data,group_new,'colors',colors_cluster,'symbol', '.')
set(gca,'FontSize',8)
xlim([0.5 4+0.5])
xlabels = {'Bit Decrease','Bit Increase','Un-recovered','Un-assigned'};
set(gca,'Xtick',1:4)
set(gca,'XtickLabel',xlabels,'FontSize',10,'FontName','Arial')
xtickangle(30)
ylabel({'Relative percent of shift','abs[(CNO-Ctrl)/Ctrl - (Pctrl-Ctrl)/Ctrl]*100'},'FontSize',10,'FontName','Arial')
ytickformat('%g%%');

p = ranksum(data(group_new == 1), data(group_new == 4))
text(0.5,3*100,['P=' num2str(p,'%.3e')],'FontSize',8)
p = ranksum(data(group_new == 2), data(group_new == 4))
text(1.5,3*100,['P=' num2str(p,'%.3e')],'FontSize',8)
p = ranksum(data(group_new == 3), data(group_new == 4))
text(2.5,3*100,['P=' num2str(p,'%.3f')],'FontSize',8)

[~,p] = kstest(data(group_new == 1))
[~,p] = ttest2(data(group_new == 1), data(group_new == 4),'Vartype','unequal')
text(0.5,3*100,['P=' num2str(p,'%.3e')],'FontSize',8)
[~,p] = ttest2(data(group_new == 2), data(group_new == 4),'Vartype','unequal')
text(1.5,3*100,['P=' num2str(p,'%.3e')],'FontSize',8)
[~,p] = ttest2(data(group_new == 3), data(group_new == 4),'Vartype','unequal')
text(2.5,3*100,['P=' num2str(p,'%.3f')],'FontSize',8)



data = abs(data_InfoScoreAll_diffPercent_relative);
data = data(idx);
shift_rel = zeros(4,length(mouse));
for i = 1:length(mouse)
    shift_rel(:,i) = grpstats(data(InfoScoreAll.mouse == i),group_new(InfoScoreAll.mouse == i),'mean');
end

data = shift_rel'*100;
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
xlabels = {'Bit Decrease','Bit Increase','Un-recovered','Un-assigned'};
set(gca,'Xtick',1:size(data,2))
set(gca,'XtickLabel',xlabels,'FontSize',12,'FontName','Arial')
xtickangle(30)
ylabel({'Relative percent of shift','abs[(CNO-Ctrl)/Ctrl - (Pctrl-Ctrl)/Ctrl]*100'},'FontSize',12,'FontName','Arial')
set(gca,'FontSize',12,'FontName','Arial')
hold on
errorbar(ctr, ydt, std(data), '.k','marker', 'none','LineWidth',1)
hold off
box off
ytickformat('%g%%');
[~,p] = ttest(shift_rel(3,:),shift_rel(4,:))
text(1.5,1.6*100,['P=' num2str(p,'%.2f')],'FontSize',8)








names_group = {'Bit Decrease','Bit Increase','Un-recovered'};
figure
for i = 1:4
    subplot(4,1,i)
    if i == 1
        cdfplot(abs(data_InfoScoreAll_diffPercent_relative)*100)
    else
        cdfplot(abs(data_InfoScoreAll_diffPercent_relative(group_new == i-1))*100)
    end
    xtickformat('%g%%');
    grid on
    set(gca,'xtick',[10 25 50 100 200 300])
    box off
    ylabel('Cumulative distribution')
    if i == 4
        xlabel({'Relative percent of shift','[(CNO-Ctrl)/Ctrl - (Pctrl-Ctrl)/Ctrl]*100'})
    else
        xlabel('')
    end
    if i == 1
        title('Distribution of all the place cells')
    else
        title(['Distribution in ',names_group{i-1}])
    end
end
data_InfoScoreAll_diffPercent_DI = data_InfoScoreAll_diffPercent;
data_InfoScoreAll_diffPercent_DI(group_new == 3,:) = [];

% determine the critical values of relative shift
Q = quantile(abs(data_InfoScoreAll_diffPercent_relative),0.05);
nboot = 1000;
data_InfoScoreAll_diffPercent_relativeboot = zeros(size(data_InfoScoreAll,1),nboot);
for i = 1:nboot
    data_InfoScoreAllboot = zeros(size(data_InfoScoreAll,1),3);
    % idx = randi(size(data_InfoScoreAll,1),1,size(data_InfoScoreAll,1));
    idx = randperm(size(data_InfoScoreAll,1));
    data_InfoScoreAllboot(:,1) = data_InfoScoreAll(idx,1);
    %  idx = randi(size(data_InfoScoreAll,1),1,size(data_InfoScoreAll,1));
    idx = randperm(size(data_InfoScoreAll,1));
    data_InfoScoreAllboot(:,2) = data_InfoScoreAll(idx,2);
    %   idx = randi(size(data_InfoScoreAll,1),1,size(data_InfoScoreAll,1));
    idx = randperm(size(data_InfoScoreAll,1));
    data_InfoScoreAllboot(:,3) = data_InfoScoreAll(idx,3);
    data_InfoScoreAll_diffPercentboot= [data_InfoScoreAllboot(:,3)./data_InfoScoreAllboot(:,1) - 1, data_InfoScoreAllboot(:,2)./data_InfoScoreAllboot(:,1) - 1];
    data_InfoScoreAll_diffPercent_relativeboot(:,i) = data_InfoScoreAll_diffPercentboot(:,2)-data_InfoScoreAll_diffPercentboot(:,1);
end
Q = quantile(abs(data_InfoScoreAll_diffPercent_relativeboot),0.95,2);
Q = quantile((data_InfoScoreAll_diffPercent_relativeboot(:)),0.95)

a = abs(data_InfoScoreAll_diffPercent_relative) < 5.57;

% overlay CNO and saline control experiments

data_InfoScoreAll_diffPercent = [data_InfoScoreAll(:,3)./data_InfoScoreAll(:,1) - 1, data_InfoScoreAll(:,2)./data_InfoScoreAll(:,1) - 1];
data_InfoScoreAll_diffPercent_relative = data_InfoScoreAll_diffPercent(:,2)-data_InfoScoreAll_diffPercent(:,1);

data_InfoScoreAll_diffPercent_relative_CNO = data_InfoScoreAll_diffPercent_relative;
group_new_CNO = group_new;
data_InfoScoreAll_diffPercent_relative_SalineCtrl = data_InfoScoreAll_diffPercent_relative;
group_new_SalineCtrl = group_new;

save data_InfoScoreAll_diffPercent.mat data_InfoScoreAll_diffPercent_relative_CNO data_InfoScoreAll_diffPercent_relative_SalineCtrl group_new_CNO group_new_SalineCtrl

names_group = {'Bit Decrease','Bit Increase','Un-recovered'};
figure
for i = 1:4
    subplot(4,1,i)
    if i == 1
        x1 = abs(data_InfoScoreAll_diffPercent_relative_CNO)*100;
        
        x2 = abs(data_InfoScoreAll_diffPercent_relative_SalineCtrl)*100;
        
    else
        x1 = abs(data_InfoScoreAll_diffPercent_relative_CNO(group_new_CNO == i-1))*100;
        x2 = abs(data_InfoScoreAll_diffPercent_relative_SalineCtrl(group_new_SalineCtrl == i-1))*100;
    end
    %         x1 = x1(x1 > 50);
    %         x2 = x2(x2 > 50);
    cdfplot(x1)
    hold on
    cdfplot(x2)
    [~,p] = kstest2(x1,x2)
    text(100,0.5,['pvalue = ',num2str(p,'%.3f')])
    %  title(['CNO:', num2str(length(x1)), ' vs saline:', num2str(length(x2))])
    xtickformat('%g%%');
    grid on
    set(gca,'xtick',[10 25 50 100 200 300])
    box off
    ylabel('Cumulative distribution')
    if i == 4
        xlabel({'Relative percent of shift','[(CNO-Ctrl)/Ctrl - (Pctrl-Ctrl)/Ctrl]*100'})
    else
        xlabel('')
    end
    if i == 1
        title(['All the place cells, CNO:', num2str(length(x1)), ' vs saline:', num2str(length(x2))] )
        legend({'CNO experiment','Saline-ctrl experiment'})
    else
        title([names_group{i-1}, ', CNO:', num2str(length(x1)), ' vs saline:', num2str(length(x2))])
    end
end


data_InfoScoreAll_diff = [data_InfoScoreAll(:,1) - data_InfoScoreAll(:,3),data_InfoScoreAll(:,1) - data_InfoScoreAll(:,2)];
data_InfoScoreAll_diffPercent = 100*[data_InfoScoreAll(:,3)./data_InfoScoreAll(:,1) - 1, data_InfoScoreAll(:,2)./data_InfoScoreAll(:,1) - 1];
figure
subplot(1,3,1)
boxplot(abs(data_InfoScoreAll_diffPercent(group_new_CNO == 1,1)))
ylabel('Percent shift: [abs(Pctrl-Ctrl)/Ctrl]*100')
xlabel('abs(Pctrl-Ctrl)')
subplot(1,3,2)
data = data_InfoScoreAll_diffPercent(group_new_CNO == 1,1);
data = data(data < 0);
boxplot(data)
xlabel('Pctrl-Ctrl < 0')
ylabel('Percent shift: [(Pctrl-Ctrl)/Ctrl]*100')
subplot(1,3,3)
data = data_InfoScoreAll_diffPercent(group_new_CNO == 1,1);
data = data(data > 0);
boxplot(data)
xlabel('Pctrl-Ctrl > 0')
ylabel('Percent shift: [(Pctrl-Ctrl)/Ctrl]*100')

figure
subplot(1,3,1)
boxplot(abs(data_InfoScoreAll_diffPercent(group_new_CNO == 1,2)))
ylabel('Percent shift: [abs(CNO-Ctrl)/Ctrl]*100')
xlabel('abs(CNO-Ctrl)')
subplot(1,3,2)
data = data_InfoScoreAll_diffPercent(group_new_CNO == 1,2);
data = data(data < 0);
boxplot(data)
xlabel('CNO-Ctrl < 0')
ylabel('Percent shift: [(CNO-Ctrl)/Ctrl]*100')
subplot(1,3,3)
data = data_InfoScoreAll_diffPercent(group_new_CNO == 1,2);
data = data(data > 0);
boxplot(data)
xlabel('CNO-Ctrl > 0')
ylabel('Percent shift: [(CNO-Ctrl)/Ctrl]*100')


% saline-ctrl
data_InfoScoreAll_diff = [data_InfoScoreAll(:,1) - data_InfoScoreAll(:,3),data_InfoScoreAll(:,1) - data_InfoScoreAll(:,2)];
data_InfoScoreAll_diffPercent = 100*[data_InfoScoreAll(:,3)./data_InfoScoreAll(:,1) - 1, data_InfoScoreAll(:,2)./data_InfoScoreAll(:,1) - 1];

figure
subplot(1,3,1)
data = abs(data_InfoScoreAll_diffPercent(group_new_SalineCtrl == 1,1));
boxplot(data)
ylabel('Percent shift: [abs(Ctrl3-Ctrl1)/Ctrl1]*100')
xlabel('abs(Ctrl3-Ctrl1)')
title(['#',num2str(nnz(group_new_SalineCtrl == 1)),'place cells'])
text(0.7,median(data),num2str(median(data),'%.0f'))
subplot(1,3,2)
data = data_InfoScoreAll_diffPercent(group_new_SalineCtrl == 1,1);
data = data(data < 0);
boxplot(data)
xlabel('Ctrl3-Ctrl1 < 0')
ylabel('Percent shift: [(Ctrl3-Ctrl1)/Ctrl1]*100')
title(['#',num2str(length(data)),'place cells'])
text(0.7,median(data),num2str(median(data),'%.0f'))
subplot(1,3,3)
data = data_InfoScoreAll_diffPercent(group_new_SalineCtrl == 1,1);
data = data(data > 0);
boxplot(data)
xlabel('Ctrl3-Ctrl1 > 0')
ylabel('Percent shift: [(Ctrl3-Ctrl1)/Ctrl1]*100')
title(['#',num2str(length(data)),'place cells'])
text(0.7,median(data),num2str(median(data),'%.0f'))

figure
subplot(1,3,1)
data = abs(data_InfoScoreAll_diffPercent(group_new_SalineCtrl == 1,2));
boxplot(data)
ylabel('Percent shift: [abs(Ctrl2-Ctrl1)/Ctrl1]*100')
xlabel('abs(Ctrl2-Ctrl1)')
title(['#',num2str(nnz(group_new_SalineCtrl == 1)),'place cells'])
text(0.7,median(data),num2str(median(data),'%.0f'))
subplot(1,3,2)
data = data_InfoScoreAll_diffPercent(group_new_SalineCtrl == 1,2);
data = data(data < 0);
boxplot(data)
xlabel('Ctrl2-Ctrl1 < 0')
ylabel('Percent shift: [(Ctrl2-Ctrl1)/Ctrl1]*100')
title(['#',num2str(length(data)),'place cells'])
text(0.7,median(data),num2str(median(data),'%.0f'))
subplot(1,3,3)
data = data_InfoScoreAll_diffPercent(group_new_SalineCtrl == 1,2);
data = data(data > 0);
boxplot(data)
xlabel('Ctrl2-Ctrl1 > 0')
ylabel('Percent shift: [(Ctrl2-Ctrl1)/Ctrl1]*100')
title(['#',num2str(length(data)),'place cells'])
text(0.7,median(data),num2str(median(data),'%.0f'))



