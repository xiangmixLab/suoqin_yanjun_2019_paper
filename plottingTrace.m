%% This function is used to overlay the neuron activity onto behaviors
% Inputs:
%        (1) neuron: a source2D variable, including identified neurons with traces and spatial information, which is obtained by runing cnmfe codes
%        (3) CellID: a vector, e.g,1:10 (display the traces of the first 10 identified neurons)
%        (3) combined: a logical variable. if it is 1 or true, it means
%        that the data consist of different conditions.

function plottingTrace(neuron,cellID,combined)
% close all
% folderName = 'FiguresCellTrace';
% if ~exist(folderName,'dir')
%     mkdir(folderName);
% end

folderName = fullfile('results','figures','FiguresCellTrace');
if ~exist(folderName, 'dir')
    mkdir(folderName);
end

fpath=folderName;
k = 0;
for tt = 1:length(cellID)
    i = cellID(tt);
    if k == 0 | mod(i-1,20) == 0
        i
        ax = figure;
        set(ax, 'Position', [100, 100, 450, 1000]);
        k = 0;
        ha = tight_subplot(10,2,[.02 .05],[.03 .02],[.05 .01]);
    end
    k = k+1;
%     subplot(5,2,k)
    axes(ha(k))
    plot(neuron.C(i,:))
%     hold on
%     plot(neuron.S(i,:),'r')
    axis tight
    if combined
        for ii = 2:length(neuron.num2read)-1
            line([sum(neuron.num2read(2:ii)) sum(neuron.num2read(2:ii))],get(gca,'Ylim'),'LineStyle','--','Color','k')
        end
    end
    set(gca,'Xtick',[]);set(gca,'FontSize',8)
    axis tight
    if mod(i,20) == 0 || i == max(cellID)
        if combined
            xtick0(1) = 1;
            for ii = 2:length(neuron.num2read)
                xtick0(ii) = sum(neuron.num2read(2:ii));
            end
            set(gca,'Xtick',xtick0,'FontSize',8)
        else
            set(gca,'Xtick',[1 ceil(neuron.num2read/2) neuron.num2read])
        end
    end
    xtickformat('%.0f')
    title(['Cell',num2str(i)],'FontName','Arial','FontSize',10)
    if mod(i,20) == 0 || i == max(cellID)
        %saveas(gcf,fullfile(fpath,strcat('TraceCell',num2str(i),'.fig')))
        saveas(gcf,fullfile(fpath,strcat('TraceCell',num2str(i),'.pdf')))
    end
end
