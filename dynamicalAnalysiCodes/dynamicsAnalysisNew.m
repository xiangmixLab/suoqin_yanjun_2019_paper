%%% temporal dynamical analysis of Ca2+ imaging data
addpath(genpath('dynamicalAnalysiCodes'))
%% load the combined neuron data file
load('1019_1021combinedNeuron.mat')
if length(neuron.num2read) == 1
    neuron.num2read = [26952, 13472, 13480]; % the first number is the total number of frames,
    % the second and third number is the number of frame in the first condition
    % and the second condition, respectively.
end
%% given a threshold for each neuron.firing
thresh = determiningFiringEventThresh(neuron); % 10% of the peak
%% filter the fake neurons for a given threshold, here is based on the peak value of trace
neuronDeleted = [];
for i = 1:size(neuron.S,1)
    if max(neuron.trace(i,:)) < 50
        neuronDeleted = [neuronDeleted,i];
    end
end
thresh(neuronDeleted) = [];
%% filter the trail if the peak is extremely small when doing clustering
threshCluster = 0.01*max(neuron.trace,[],2); % a vector
threshCluster(neuronDeleted) = [];
%% load the individual neuron data and extract the time information in each condition
load('neuronIndividuals.mat')
downsampeFactor = round(length(neuronIndividuals{1}.time)/size(neuronIndividuals{1}.trace,2));
if downsampeFactor == 2
    for i = 1:length(neuronIndividuals)
        neuronIndividuals{i}.time = neuronIndividuals{i}.time(1:downsampeFactor:end);
    end
else
    disp('Please check the downsampeFactor')
end
% delete the fake neurons
for i = 1:length(neuronIndividuals)
    neuronIndividuals{i}.S(neuronDeleted,:) = [];
    neuronIndividuals{i}.trace(neuronDeleted,:) = [];
    neuronIndividuals{i}.centroid(neuronDeleted,:) = [];
    neuronIndividuals{i}.C(neuronDeleted,:) = [];
    neuronIndividuals{i}.A(:,neuronDeleted) = [];
    neuronIndividuals{i}.Coor(neuronDeleted) = [];
end
neuron1 = neuronIndividuals{1};
neuron2 = neuronIndividuals{2};

%% extract each trail
% update the neuron information;
% idxLRi is the trail information from left to right, each column is one trail
% idxRLi is the trail information from right to right, each column is one trail
% idxTrail1 is all the identified trail
behav1 = importdata('BehavCA1_1019linearNoCNO.mat');
[neuron1,idxLRi1,idxRLi1,idxTrail1] = extractTrailLinear(neuron1, behav1);

behav2 = importdata('BehavCA1_1021_linearCNO.mat');
[neuron2,idxLRi2,idxRLi2,idxTrail2] = extractTrailLinear(neuron2, behav1); % update the neuron information

%% select the condition to be analyzed
neuron0 = neuron1;
idxLRi = idxLRi1;
idxRLi = idxRLi1;
idxTrail = idxTrail1;
%% display the trail and overlay the neuron activity onto trajectories
segDisplay = 1:6; % only dislay the first six trails in each direction
displayIndividualTrails(neuron0,segDisplay,thresh)

%% performing clustering using kmeans+consensus clustering method
% using neuron.trace
K = 5; % initial guess of the number of clusters
N = 1000; % the number of repeated times
CM = consensusKmeans(idxLRi,idxRLi,neuron0,threshCluster,K,N); % return the simimarity matrix between paired neurons
% Note: please determine the number of optimal clusters based on the pop-up clustergram
optimalK = 6; % the optimal number of clusters
Z = linkage(CM,'complete');
group = cluster(Z,'maxclust',optimalK);
threshC = Z(end-optimalK+2,3)-eps;

% reproduce the clustergram and assign different colors to each cluster
cgo = clustergram(CM,'Standardize',3,'Linkage','complete','Colormap',redbluecmap);
set(cgo,'Dendrogram',threshC);
% colorClustersC = cell(length(optimalK),1);
% for i = 1:optimalK
%     colorClustersC{i} = colorClusters(i,:);
% end
% rm = struct('GroupNumber',{2051,2057,2058,2049,2059,2056,2055,2043},'Annotation',{'I','II','III','IV','V','VI','VII','VIII'},...
%      'Color',colorClustersC);
% set(cgo,'RowGroupMarker',rm)

%% rename the cluster asignment variables which give the cluster assignment for each neuron
CMCT = CM;
groupCT = group;

CMCNO = CM;
groupCNO = group;

%%  extract the order of neuron in the clustergram
cgolabels = cgo.RowLabels;
[cgolabels,~,perm] = intersect(cgolabels,cellstr(num2str([1:length(cgolabels)]')),'stable');

permCT = perm; % please rename
permCNO = perm;

%% generate different colors
colorClusters = distinguishable_colors(2*optimalK+1);
colorClusters(4,:) = []; % the fourth is black
colorClusters(end-optimalK+1:end,:) = []; % colors for condition 1
colorClusters(1:optimalK,:) = []; % colors for condition 2

%% display the simimarity matrix using heatmap
% displaySimilarityHeatmap(CM1,CM2,perm,colorClusters)
% CM1 is for showing the desired similarity matrix;
% CM2,perm,colorClusters are for reordering. If one wants to show the
% heatmap of CNO, but reorder the neuron according to the order in Control,
% then CM1 is the similarity of CNO, CM2,perm,colorClusters are the
% variables of Control.
displaySimilarityHeatmap(CM,CM,perm,colorClusters)

%% display the trace of each cluster
trialNum = [1:10]; % the trails for displaying, e.g. the first ten trails
displayTrace(neuron0,trialNum, group,idxTrail)

%% display spatial map
% colorCell6 = [248 118 109; 183 159 0; 0 186 56;0 191 196; 97 156 255; 245 100 227]/255;
showCenter = 1;
displayspatialMap(neuron0,group,colorClusters,showCenter)


%% compare different conditions
%% calculating the percentage of neurons changed in each cluster
count = zeros(optimalK,optimalK);
for i = 1:optimalK
    count(i,:) = hist(groupCNO(groupCT == i),1:optimalK);
end
%     percent = count./sum(count,2);
figure
for i = 1:optimalK
    subplot(2,3,i)
    h = pie(count(i,:));
    idx = find(count(i,:));
    for ii = 1:length(idx)
        h(2*ii-1).FaceColor = colorClusters(idx(ii),:);
    end
    title(['n=',num2str(sum(count(i,:)))]);
    if i == 2
        labels = cellstr(num2str([1:optimalK]'));
        legend(labels,'Location','southoutside','Orientation','horizontal')
    end
end

%% detect overlap clusters in two conditions
neuron22trace= zeros(size(neuron2.trace,1),length(neuron1.time));
for i = 1:size(neuron2.trace,1)
    neuron22trace(i,:) = interp1(neuron2.time, neuron2.trace(i,:)', neuron1.time,'linear','extrap');
end
D = zeros(length(unique(groupCT)),length(unique(groupCNO)));
for i = 1:length(unique(groupCT))
    for j = 1:length(unique(groupCNO))
        %         d = pdist2(neuron1.trace(groupCT == i,:),neuron22trace(groupCNO == j,:),'correlation');
        %         d = corr(neuron1.trace(groupCT == i,:)',neuron22trace(groupCNO == j,:)');
        d = corr(mean(neuron1.trace(groupCT == i,:))',mean(neuron22trace(groupCNO == j,:))');
        D(i,j) = mean(d(:));
    end
end
% display the correlation between pairwise clusters in each condition
labels = cellstr(num2str([1:optimalK]'));
figure;
heatmap1(D, labels, labels,'%0.2f', 'Colormap', 'money',...
    'Colorbar', true, 'FontSize', 8,'ShowAllTicks', true);
% c = max(abs([min(PCCTF_TF(:)),max(PCCTF_TF(:))]));
% caxis([-c c])
c = colorbar;
c.Location = 'east';
 
c.FontSize = 8;
xlabel('CNO clusters','FontSize',10);ylabel('Control clusters','FontSize',10)
title('Correlations of pairwise clusters','FontSize',10)


% % neuron22S= zeros(size(neuron2.S,1),length(neuron1.time));
% % for i = 1:size(neuron2.S,1)
% % neuron22S(i,:) = interp1(neuron2.time, neuron2.S(i,:)', neuron1.time,'linear','extrap');
% % end
% % DS = zeros(length(unique(groupCT)),length(unique(groupCNO)));
% % for i = 1:length(unique(groupCT))
% %     for j = 1:length(unique(groupCNO))
% % %         d = pdist2(neuron1.trace(groupCT == i,:),neuron22trace(groupCNO == j,:),'correlation');
% % %         d = corr(neuron1.trace(groupCT == i,:)',neuron22trace(groupCNO == j,:)');
% %         d = corr(mean(neuron1.S(groupCT == i,:))',mean(neuron22S(groupCNO == j,:))');
% %         DS(i,j) = mean(d(:));
% %     end
% % end

%% show the dynamics of each cluster

data1C = [];data2C = [];
for i = 1:length(unique(groupCT))
    data1C = [data1C;neuron1.S(groupCT == i,:)];
    data2C = [data2C;neuron2.S(groupCT == i,:)];
end
position = 0;
for i = 1:optimalK
    position(i+1) = position(i)+sum(groupCT == i);
end

figure
subplot(1,2,1); imagesc(data1C)
hold on
flag = 1;
if flag
    for i = 2:length(position)-1
        line(get(gca,'XLim'),[position(i) position(i)],'LineWidth',0.5,'Color','k'); hold on;
        hold on;
    end
end
c = colorbar;
% c.Location = 'northoutside';
% c.Label.String = 'Trace';
% c.Label.FontSize = 8;c.Label.FontWeight = 'bold';c.FontSize = 6;
caxis([0 max([data1C(:);data2C(:)])])
set(gca,'FontSize',8)
labels = strcat('C',cellstr(num2str([1:optimalK]')));
xtick0 = position(1:end-1)+diff(position)/2;
set(gca,'Ytick',xtick0);set(gca,'YtickLabel',labels,'FontName','Arial','FontSize',10)
xlabel('Frames','FontSize',10)
colormap(flipud(hot))

subplot(1,2,2); imagesc(data2C)
hold on
flag = 1;
if flag
    for i = 2:length(position)-1
        line(get(gca,'XLim'),[position(i) position(i)],'LineWidth',0.5,'Color','k'); hold on;
        hold on;
    end
end
c = colorbar;
% c.Location = 'northoutside';
% c.Label.String = 'Trace';
% c.Label.FontSize = 8;c.Label.FontWeight = 'bold';c.FontSize = 6;
caxis([0 max([data1C(:);data2C(:)])])
c.Ticks = [];
set(gca,'FontSize',8)
% labels = strcat('C',cellstr(num2str([1:optimalK]')));
% xtick0 = position(1:end-1)+diff(position)/2;
% set(gca,'Ytick',xtick0);set(gca,'YtickLabel',labels,'FontName','Arial','FontSize',8)
set(gca,'Ytick',[])
xlabel('Frames','FontSize',10)
colormap(flipud(hot))

suptitle('Trace patterns');
%% intra-cluster and inter-cluster pairwise cell distance analysis
% (1) intra-cluster distance analysis
% (1.1) comparison the distance between Control and CNO conditions
DE = zeros(length(unique(groupCT)),2);
DC = zeros(length(unique(groupCT)),2);
for i = 1:length(unique(groupCT))
    d = pdist(neuron1.trace(groupCNO == i,:));
    DE(i,1) = mean(d(:));
    d = pdist(neuron2.trace(groupCNO == i,:));
    DE(i,2) = mean(d(:));
    
    d = pdist(neuron1.trace(groupCNO == i,:),'correlation');
    DC(i,1) = mean(1-d(:));
    d = pdist(neuron2.trace(groupCNO == i,:),'correlation');
    DC(i,2) = mean(1-d(:));
end
% display the average distance across all the cells in an individual cluster
figure
subplot(1,2,1);bar(DE)
legend({'Control','CNO'});
title('Cells grouped by Control clusters','FontSize',10,'FontWeight','bold')
xlabel('Clusters','FontSize',10);
ylabel('Avg intra-cluster pairwise cell distance','FontSize',10)
xlim([0.5 size(DE,1)+0.5])
subplot(1,2,2);bar(DC)
legend({'Control','CNO'});
title('Cells grouped by Control clusters','FontSize',10,'FontWeight','bold')
xlabel('Clusters','FontSize',10);
ylabel('Avg intra-cluster pairwise correlations','FontSize',10)
xlim([0.5 size(DC,1)+0.5])

% display the average distance across all the clusters
figure
subplot(1,2,1);bar(mean(DE))
set(gca,'XtickLabel',{'Control','CNO'})
ylabel('Avg intra-cluster pairwise cell distance','FontSize',10)
xlim([0.5 2.5])
subplot(1,2,2);bar(mean(DC))
set(gca,'XtickLabel',{'Control','CNO'})
ylabel('Avg intra-cluster pairwise correlations','FontSize',10)
xlim([0.5 2.5])
suptitle('Cells grouped by Control clusters');

