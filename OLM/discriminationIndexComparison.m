function discriminationIndexComparison(DIall_saline,DIall_CNO,sessionID,sessionName)
folderName = fullfile('results','figures');
if ~exist(folderName, 'dir')
    mkdir(folderName);
end
measureNames = {'sumFiringRate','sumCountTime','sumCount','sumAmplitude'};
measures = {'Firing rate','Exploration time','Firing count','Amplitude'};
conditions = {'Saline','CNO'};
colors_conditions = [228,26,28;55,126,184]/255;
hFig = figure('position', [200, 200, 680,150]);
% ha = tight_subplot(1,3,[.001 .1],[.1 .1],[.1 .02]);
for j  = 1:length(measures)
    data1 = DIall_saline.(char(measureNames(j)))(:,sessionID);
    data2 = DIall_CNO.(char(measureNames(j)))(:,sessionID);
    group = ones(size(data1,1)+size(data2,1),1)*2;group(1:size(data1,1)) = 1;
    data = [data1(:);data2(:)];
    data(isnan(data)) = 0;
    hF = subplot(1,4,j);
    %     axes(ha(j));
    %     boxplot(data,group,'symbol', ' ','colors',colors_session)
    boxplot(data,group,'symbol', ' ')
    hold on
    plotSpread({data1,data2},'distributionMarkers', {'o', 'o'},'distributionColors','k')
    xlim([0.5 2+0.5])
    set(gca,'Xtick',1:2)
    set(gca,'XtickLabel',conditions,'FontSize',12,'FontName','Arial')
    ylabel('Discrimination index','FontSize',12,'FontName','Arial')
    title(measures{j},'FontSize',12,'FontName','Arial')
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
    %   suptitle(sessionName)
    
    show_pvalue = 1;
    if show_pvalue
        [~,p1] = kstest(data(group == 1))
        [~,p2] = kstest(data(group == 2))
        %         [~,p] = ranksum(data(group == 1),data(group == 2),'tail','right')
        if p1 < 0.05 & p2 < 0.05
%             if strcmpi(sessionName,'testing')
%                 [~,p] = ttest2(data(group == 1),data(group == 2),'tail','right');
%             elseif strcmpi(sessionName,'training')
%                 [~,p] = ttest2(data(group == 1),data(group == 2));
%             else
%                 sprintf('Please input an correct session name: training or testing')
%             end
            [~,p] = ttest2(data(group == 1),data(group == 2));
        else
            p = anova1(data,group,'off');
        end
        
        %  p = anova1(data,group,'off');
        pos = get(gca,'position');dim = [pos(1)+pos(3)/4,max(0,pos(2)-pos(4)/30),pos(3),pos(4)];
        
        if p < 0.01
            annotation('textbox',dim,'String',['P = ', num2str(p,'%3.0e')],'FitBoxToText','on','LineStyle','none','EdgeColor','none','FontSize',8,'FontWeight','bold');
        else
            annotation('textbox',dim,'String',['P = ', num2str(p,'%.3f')],'FitBoxToText','on','LineStyle','none','EdgeColor','none','FontSize',8,'FontWeight','bold');
        end
    end
    ylim([min(data)-0.05 max(data)+0.1])
    saveas(gcf,fullfile(folderName,['discriminationIndexComparison',sessionName,'.pdf']))
    saveas(gcf,fullfile(folderName,['discriminationIndexComparison',sessionName,'.fig']))
end