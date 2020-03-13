function displayspatialMap(neuron0,group,colorClusters,showCenter,showShape)
colorCell = zeros(length(group),3);
for i = 1:length(group)
    colorCell(i,:) = colorClusters(group(i),:);
end
figure
% showCenter = 1;
if showCenter
    %gscatter(neuron0.centroid(:,2),neuron0.centroid(:,1),group,colorClusters,[],6)
    for i = 1:length(unique(group))
        h(i) = scatter(neuron0.centroid(group == i,2),neuron0.centroid(group == i,1),8,colorClusters(i,:),'filled');
        hold on
    end
    box on
    
%     legend('C1','C2','C3','C4')
end
hold on
if showShape
Coor1 = neuron0.Coor;
for i = 1:length(Coor1)
    cont = medfilt1(Coor1{i}')';
    if size(cont,2) > 1
        plot(cont(1,2:end),cont(2,2:end),'Color',colorCell(i,:), 'linewidth', 1.5); hold on;
    end
end
hold on
end
axis ij
set(gca,'color','w')
% set(gca,'Xtick',[]);set(gca,'Ytick',[]);
% axis([0 max(neuron0.pos(:,1)) 0 max(neuron0.pos(:,2))])

