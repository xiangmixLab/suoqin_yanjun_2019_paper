%%% temporal dynamical analysis of Ca2+ imaging data
addpath(genpath('dynamicalAnalysiCodes'))
%% load the combined neuron data file
addpath(genpath('E:\Google Drive synchronization\projects\Project2_Xu\Miniscope_Matlab_code'))
filefolder = 'H:\\Miniscope data analysis_Yanjun\\M3244F';
%% load the combined neuron data file
% load('1019_1021combinedNeuron.mat')
load(fullfile(filefolder,'M3244F_Ctrl2_CNO1_CNO2_PCtrl1_PCtrl2_Neuron2.mat'))
%% given a threshold for each neuron.firing
thresh = determiningFiringEventThresh(neuron,'S'); % 10% of the peak
%% load the individual neuron data and extract the time information in each condition
load(fullfile(filefolder,'neuronIndividuals.mat'))
%% filter the fake neurons for a given threshold, here is based on the peak value of trace
% neuronDeleted = [];
% for i = 1:size(neuron.S,1)
%     if max(neuron.trace(i,:)) < 50
%         neuronDeleted = [neuronDeleted,i];
%     end
% end
load('M3244F_0217Ctrl_InfoPerSecAll.mat')
idx_PC = find(infoPerSecond < 0.15);
neuronDeleted = idx_PC;
thresh(neuronDeleted) = [];
%% filter the trail if the peak is extremely small when doing clustering
threshCluster = 0.01*max(neuron.trace,[],2); % a vector
threshCluster(neuronDeleted) = [];
%% load the individual neuron data and extract the time information in each condition
% load('neuronIndividuals.mat')
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
mice = 'M3244F';section = {'Contrl2','CNO1','CNO2','Post-contrl1','Post-contrl2'};
behavIndividuals = cell(1,length(neuronIndividuals));
for i = 1:length(neuronIndividuals)
    behavIndividuals{i} = importdata([mice,'_',section{i},'_Behav.mat']);
end

idxLRIndi = cell(1,length(neuronIndividuals));idxRLIndi = idxLRIndi;idxTrailIndi = idxRLIndi;
for i = 1:length(neuronIndividuals)
    [neuron2,idxLRi2,idxRLi2,idxTrail2] = extractTrailLinear(neuronIndividuals{i}, behavIndividuals{i});
    neuronIndividuals{i} = neuron2;
    idxLRIndi{i} = idxLRi2;idxRLIndi{i} = idxRLi2;idxTrailIndi{i} = idxTrail2;
end
%% select the condition to be analyzed
neuronID = 1;
neuron0 = neuronIndividuals{neuronID};
idxLRi = idxLRIndi{neuronID};
idxRLi = idxRLIndi{neuronID};
idxTrail = idxTrailIndi{neuronID};
%% display the trail and overlay the neuron activity onto trajectories
segDisplay = 1:6; % only dislay the first six trails in each direction
displayIndividualTrails(neuron0,segDisplay,thresh)

%% performing clustering using kmeans+consensus clustering method
% using neuron.trace
K = 3; % initial guess of the number of clusters
N = 50; % the number of repeated times
%CM = consensusKmeans(idxLRi,idxRLi,neuron0,threshCluster,K,N); % return the simimarity matrix between paired neurons
CM = consensusKmeans(neuron0,threshCluster,K,N); 
% Note: please determine the number of optimal clusters based on the pop-up clustergram
optimalK = 3; % the optimal number of clusters
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

CMPC1 = CM;
groupPC1 = group;

%%  extract the order of neuron in the clustergram
cgolabels = cgo.RowLabels;
[cgolabels,~,perm] = intersect(cgolabels,cellstr(num2str([1:length(cgolabels)]')),'stable');

permCT = perm; % please rename
permCNO = perm;
permPC1 = perm;

%% generate different colors
colorClusters = distinguishable_colors(optimalK+1);
colorClusters(4,:) = []; % the fourth is black
colorClusters  = cmap(1:3,:);
%% display the simimarity matrix using heatmap
% displaySimilarityHeatmap(CM1,CM2,perm,colorClusters)
% CM1 is for showing the desired similarity matrix;
% CM2,perm,colorClusters are for reordering. If one wants to show the
% heatmap of CNO, but reorder the neuron according to the order in Control,
% then CM1 is the similarity of CNO, CM2,perm,colorClusters are the
% variables of Control.
displaySimilarityHeatmap(CM,CM,perm,colorClusters)
displaySimilarityHeatmap(CMPC1,CMCT,permCT,colorClusters) % order neurons in the heatmap of post-control (PC) case based on the order of control (permCT) 
%% display the trace of each cluster
trialNum = [1:10]; % the trails for displaying, e.g. the first ten trails
displayTrace(neuron0,trialNum, groupCNO,idxTrail,colorClusters)

%% display spatial map
% colorCell6 = [248 118 109; 183 159 0; 0 186 56;0 191 196; 97 156 255; 245 100 227]/255;
showCenter = 1;showShape = 0;
displayspatialMap(neuron0,group,colorClusters,showCenter,showShape)
axis([0 250 0 160])
set(gca,'Xtick',[]);set(gca,'Ytick',[]);
save variables_M3244F_Ctrl2_CNO1_CNO2_PCtrl1_PCtrl2_clustering_0.35PC.mat

%% using heatmap to show the firing/trace dynamics with the changes of frames/time
neuron.S(neuronDeleted,:) = [];
neuron.trace(neuronDeleted,:) = [];
dataC = [];
for j = 1:length(unique(groupCT))
    dataC = [dataC;neuron.trace(groupCT == j,:)];
end

positionS = 0;
for i = 1:length(neuronIndividuals)
    positionS(i+1) = positionS(i)+neuron.num2read(i+1);
    dataCi = [];
end
positionC = 0;
for i = 1:optimalK
    positionC(i+1) = positionC(i)+sum(groupCT == i);
end

figure
imagesc(dataC)
colormap(flipud(hot))
colorbar
hold on
flag = 1;
if flag
    for i = 2:length(positionS)-1
        line([positionS(i)+0.5 positionS(i)+0.5],get(gca,'YLim'),'LineWidth',0.5,'Color','k'); hold on;
    end
    for i = 2:length(positionC)-1
        line(get(gca,'XLim'),[positionC(i) positionC(i)],'LineWidth',0.5,'Color','k'); hold on;
    end
end
set(gca,'FontSize',8)
labels = strcat('C',cellstr(num2str([1:optimalK]')));
ytick0 = positionC(1:end-1)+diff(positionC)/2;
set(gca,'Ytick',ytick0);set(gca,'YtickLabel',labels,'FontName','Arial','FontSize',10)
xtick0 = positionS(1:end-1)+diff(positionS)/2;
xlabels = {'Contrl (Day1)','CNO (Day 6)', 'CNO (Day 8)','Post-contrl (Day 10)','Post-contrl (Day 12)'};
set(gca,'Xtick',xtick0);set(gca,'XtickLabel',xlabels,'FontName','Arial','FontSize',10)
xtickangle(45)
xlabel('Frames','FontSize',10)
ylabel('Neurons','FontSize',10)

%% comparing the average trace (neuron.trace) among different conditions
dataS = zeros(size(neuronIndividuals{1}.S,1),2);
position = 0;
for i = 1:length(neuronIndividuals)
    dataS(:,i) = sum(neuronIndividuals{i}.trace,2)/neuronIndividuals{i}.num2read;
    position(i+1) = position(i)+1;
end
dataSC = [];
for j = 1:length(unique(groupCT))
    dataSC = [dataSC;dataS(groupCT == j,:)];
end
figure
imagesc(dataSC)
colormap(flipud(hot))
c = colorbar;
c.Location = 'southoutside';
c.Label.String = 'dF/F per frame';
c.Label.FontSize = 8;c.Label.FontWeight = 'bold';c.FontSize = 6;
hold on
flag = 1;
if flag
    for i = 2:length(position)-1
        line([position(i)+0.5 position(i)+0.5],get(gca,'YLim'),'LineWidth',0.5,'Color','k'); hold on
    end
    for i = 2:length(positionC)-1
        line(get(gca,'XLim'),[positionC(i) positionC(i)],'LineWidth',0.5,'Color','k'); hold on;
    end
end
set(gca,'FontSize',8)
labels = strcat('C',cellstr(num2str([1:optimalK]')));
ytick0 = positionC(1:end-1)+diff(positionC)/2;
set(gca,'Ytick',ytick0);set(gca,'YtickLabel',labels,'FontName','Arial','FontSize',10)
set(gca,'Xtick',1:length(neuronIndividuals));
xlabels = {'Contrl (Day1)','CNO (Day 6)', 'CNO (Day 8)','Post-contrl (Day 10)','Post-contrl (Day 12)'};
set(gca,'XtickLabel',xlabels,'FontName','Arial','FontSize',10)
xtickangle(45)
ylabel('Neurons','FontSize',10)

%% intra-cluster and inter-cluster pairwise cell distance analysis
% (1) intra-cluster distance analysis
% (1.1) comparison the distance between Control and CNO conditions
DE = zeros(length(unique(groupCT)),2);
DC = zeros(length(unique(groupCT)),length(neuronIndividuals));
for i = 1:length(unique(groupCT))
    %     d = pdist(neuronIndividuals{1}.trace(groupCT == i,:));
    %     DE(i,1) = mean(d(:));
    %     d = pdist(neuronIndividuals{5}.trace(groupCT == i,:));
    %     DE(i,2) = mean(d(:));
    for k = 1:length(neuronIndividuals)
        d = pdist(neuronIndividuals{k}.trace(groupCT == i,:),'correlation');
        DC(i,k) = mean(1-d(:));
    end
end
% display the average distance across all the cells in an individual cluster
% figure
% subplot(1,2,1);bar(DE)
% legend({'Control','CNO'});
% title('Cells grouped by Control clusters','FontSize',10,'FontWeight','bold')
% xlabel('Clusters','FontSize',10);
% ylabel('Avg intra-cluster pairwise cell distance','FontSize',10)
% xlim([0.5 size(DE,1)+0.5])
figure
% subplot(1,2,2);
b = bar(DC,'FaceColor','flat');
% legend({'Control','CNO'});
% title('Cells grouped by Control clusters','FontSize',10,'FontWeight','bold')
load colormap7.mat
for k = 1:size(DC,2)
    b(k).FaceColor = cmap(k,:);
end
set(gca,'XtickLabel',{'C1','C2','C3'},'FontName','Arial','FontSize',10);
xlabels = {'Ctrl','CNO','PCtrl'};
legend(xlabels,'FontName','Arial','FontSize',8)
xlabel('Clusters','FontSize',10);
ylabel('Avg intra-cluster pairwise correlations','FontSize',10)
xlim([0.5 size(DC,1)+0.5])

