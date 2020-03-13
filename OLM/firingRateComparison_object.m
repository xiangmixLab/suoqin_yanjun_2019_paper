function firingRateComparison_object(sumFiringRateObjectAllAll_saline,sumFiringRateObjectAllAll_CNO)
folderName = fullfile('results','figures');
if ~exist(folderName, 'dir')
    mkdir(folderName);
end

data1 = [];
for i = 1:size(sumFiringRateObjectAllAll_saline,1)
    data1 = [data1;sumFiringRateObjectAllAll_saline{i,2}];
end
data2 = [];
for i = 1:size(sumFiringRateObjectAllAll_CNO,1)
    data2 = [data2;sumFiringRateObjectAllAll_CNO{i,2}];
end

data3 = [];
for i = 1:size(sumFiringRateObjectAllAll_saline,1)
    data3 = [data3;sumFiringRateObjectAllAll_saline{i,3}];
end

data4 = [];
for i = 1:size(sumFiringRateObjectAllAll_CNO,1)
    data4 = [data4;sumFiringRateObjectAllAll_CNO{i,3}];
end
dataAll = {data1,data2,data3,data4};




conditions = {'Saline (Training)','CNO (Training)','Saline (Testing)','CNO (Testing)'};
xlabels = {'Unmoved object','Moved object'};
colors_conditions = [228,26,28;55,126,184]/255;
hFig = figure('position', [200, 200, 680,200]);
% ha = tight_subplot(1,3,[.001 .1],[.1 .1],[.1 .02]);
for j  = 1:length(dataAll)
    if j == 1 | j == 3
        colors_conditions = [228,26,28;228,26,28]/255;
    else
        colors_conditions = [55,126,184;55,126,184]/255;
    end
    data = dataAll{j};
    hF = subplot(1,4,j);
    %     axes(ha(j));
    %     boxplot(data,group,'symbol', ' ','colors',colors_session)
    %     boxplot(data,'symbol', ' ')
    %     hold on
    %     plotSpread({data(:,1),data(:,2)},'distributionMarkers', {'o', 'o'},'distributionColors','k')
    %     xlim([0.5 2+0.5])
    %     set(gca,'Xtick',1:2)
    %     set(gca,'XtickLabel',xlabels,'FontSize',12,'FontName','Arial')
    %     ylabel('Firing rate','FontSize',12,'FontName','Arial')
    %     title(conditions{j},'FontSize',12,'FontName','Arial')
    %     set(gca,'FontSize',12,'FontName','Arial')
    %     xtickangle(45)
    %     % Alter linestyle
    %     idxColor = fliplr(1:2);
    %     h2 = findobj(gca,'Tag','Box');
    %     for jj=1:length(h2)
    %         patch(get(h2(jj),'XData'),get(h2(jj),'YData'),colors_conditions(idxColor(jj),:),'FaceAlpha',0.6);
    %     end
    %     h3 = findobj(gca,'tag','Outliers');
    %     for jj = 1:length(h3)
    %         h3(jj).MarkerEdgeColor = colors_conditions(idxColor(jj),:);
    %     end
    %
    %     box off
    %
    %     lines = findobj(gcf,'Type','Line');
    %     set(lines,'LineWidth',1)
    %     h1 = findobj(gcf,'tag','Median');
    %     set(h1,{'linew'},{2.5})
    %     set(h1,'Color',[0 0 0])
    %     h = findall(gcf,'marker','.');
    %     set(h, 'markersize', 10)
    %
    %     set(gca,'LineWidth',1)
    %     set(gca,'FontName','Arial')
    %     set(gca,'FontSize',12)
    %
    %     txt_obj = findall(gca,'Type','text');
    %     set(txt_obj,'FontName','Arial','FontSize',12);
    %     %   suptitle(sessionName)
    
    
    hBar = bar(mean(data),0.4);
    ctr = [];ydt = [];
    for k1 = 1:length(hBar)
        ctr(k1,:) = bsxfun(@plus, hBar(1).XData, [hBar(k1).XOffset]');
        ydt(k1,:) = hBar(k1).YData;
    end
    hBar.FaceColor = 'flat';
    for jj = 1:size(hBar.CData,1)
        hBar.CData(jj,:) = colors_conditions(jj,:);
    end
    xlim([0.5 size(hBar.CData,1)+0.5])
    plotSpread({data(:,1),data(:,2)},'distributionMarkers', {'o', 'o'},'distributionColors','k')
    
    xlim([0.5 size(data,2)+0.5])
    set(gca,'Xtick',1:size(data,2))
    set(gca,'XtickLabel',xlabels,'FontSize',12,'FontName','Arial')
    xtickangle(45)
    ylabel('Ensemble event rate','FontSize',12,'FontName','Arial')
    set(gca,'FontSize',12,'FontName','Arial')
    hold on
    errorbar(ctr, ydt, std(data), '.k','marker', 'none','LineWidth',1)
    hold off
    box off
    
    title(conditions{j},'FontSize',12,'FontName','Arial')
    lines = findobj(gcf,'Type','Line');
    set(lines,'LineWidth',1)
    
    set(gca,'LineWidth',1)
    set(gca,'FontName','Arial')
    set(gca,'FontSize',12)
    
    txt_obj = findall(gca,'Type','text');
    set(txt_obj,'FontName','Arial','FontSize',12);
    
    
    show_pvalue = 1;
    if show_pvalue
        %         [~,p1] = kstest(data(:,1))
        %         [~,p2] = kstest(data(:,2))
          [~,p] = ttest(data(:,1),data(:,2));
%         if j == 1 | j == 2
%             [~,p] = ttest(data(:,1),data(:,2));
%         else
%             if mean(data(:,1)) > mean(data(:,2))
%                 [~,p] = ttest(data(:,1),data(:,2),'tail','right');
%             else
%                 [~,p] = ttest(data(:,1),data(:,2),'tail','left');
%             end
%         end
        %  p = anova1(data,group,'off');
        pos = get(gca,'position');dim = [pos(1)+pos(3)/4,max(0,pos(2)-pos(4)/30),pos(3),pos(4)];
        
        if p < 0.01
            annotation('textbox',dim,'String',['P = ', num2str(p,'%3.0e')],'FitBoxToText','on','LineStyle','none','EdgeColor','none','FontSize',8,'FontWeight','bold');
        else
            annotation('textbox',dim,'String',['P = ', num2str(p,'%.3f')],'FitBoxToText','on','LineStyle','none','EdgeColor','none','FontSize',8,'FontWeight','bold');
        end
    end
    %     ylim([min(data(:))-0.05 max(data(:))+0.1])
    ylim([0 max(data(:))+8])
    saveas(gcf,fullfile(folderName,['ensemblefiringRateComparison_object_bar.pdf']))
    saveas(gcf,fullfile(folderName,['ensemblefiringRateComparison_object_bar.fig']))
end