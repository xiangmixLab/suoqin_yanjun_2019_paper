function displaySimilarityHeatmap(CM1,CM2,perm,colorClusters)
figure
subplot(1,10,[3 8])
CM1(logical(eye(size(CM1)))) = 1;
imagesc(flip(CM1(perm,perm)));
axis square
set(gca,'Xtick',[]);set(gca,'Ytick',[])
colormap(bluewhitered2)
caxis([0 1])
c = colorbar;
c.Location = 'east';
c.Label.String = 'Simimarity';
c.Label.FontSize = 8;%c.Label.FontWeight = 'bold'
c.FontSize = 8;

subplot(1,10,[1 2])
optimalK = length(colorClusters);
Z = linkage(CM2,'complete');
group = cluster(Z,'maxclust',optimalK);

H = dendrogram(Z,0,'Orientation','left','Labels',[],'reorder',perm);  %[H,T,PERM] = dendrogram(Z,0,'colorthreshold',0.35); %dendrogram plot with colorthreshold set to 0.35
set(gca,'fontsize',8);
set(H,'color','k')
ylim([0.5 length(group)+0.5])
set(gca,'Xtick',[]);set(gca,'Ytick',[]);

subplot(1,10,[9.7 10])
groupUni = unique(group(perm),'stable');
ybar = [];
for i = 1:1:optimalK
    ybar = [ybar,nnz(group == groupUni(i))];
end
H = bar([ybar;nan(1,length(ybar))],'stacked','BarWidth',0.0004,'ShowBaseLine','off','LineWidth',0.01);
k = 0;
for i = 1:1:optimalK
    k = k+1;
    set(H(i),'facecolor',colorClusters(groupUni(k),:))
end
axis tight
xlim([1 1.001])
set(gca,'xtick',1)
axis off
