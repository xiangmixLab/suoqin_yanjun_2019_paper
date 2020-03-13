function [firingRateAll,countAll,countTime] = calculatingCellSpatialLinearTrackForSingleData2(neuron,behav,Segment,threshold,temp,intensity,plotting)
%% Inputs:
%     neuron: a source2D variable, including identified neurons with traces and spatial information, which is obtained by runing cnmfe codes
%     behav: behavior information, which is obtained by using Tristan's code
% Segment: a vector, e.g,1:10 (display the traces of the first 10 identified neurons)
% threshold:the threshold above which the neuron is active, e.g.,0.1
% temp:judge whether the neuron is active using neuron.trace or neuron.S; temp = 'trace' ot temp = 'S'
% downsampling: if downsampling = true, then do downsampling for neuron.time; default is false
% intensity: if intensity = true, then diplay the peak map; otherwise display the heatmap of firating rate; default is false
%%% example usage: e.g.1. plottingCellSpatialForSingleData(neuron,behav,1:5)
%%%                e.g.2, plottingCellSpatialForSingleData(neuron,behav,1:5,0.1,'trace',true,true)
%%%                e.g.3, plottingCellSpatialForSingleData(neuron,behav,1:5,0.1,'trace')
if nargin<7;     plotting = false; end
if nargin<6;     intensity = false; end
if nargin<5;     temp = 'S'; end
if nargin<4;     threshold = 0.1; end
if nargin<3;     Segment = 1:size(neuron.trace,1); end

downsampling = length(neuron.time)/size(neuron.trace,2);
if downsampling ~= 1
    %     downsampling == 2
    neuron.time = double(neuron.time);
    neuron.time = neuron.time(1:downsampling:end);
end
t = find(diff(behav.time)<=0);
behav.time(t+1) = behav.time(t)+1;
neuron.pos = interp1(behav.time,behav.position,neuron.time); %%


folderName = 'FiguresCellSpatial';
if ~exist(folderName,'dir')
    mkdir(folderName);
end
fpath=folderName;
global ts pos1 pos2

% pos1 = 0:10:ceil(max(neuron.pos(:,1)));pos2 = 0:10:ceil(max(neuron.pos(:,2)));
% pos1 = 0:binsize:ceil(max(neuron.pos(:,1)));pos2 = 0:binsize:ceil(max(neuron.pos(:,2)));

num = length(Segment);
edge11 = max(neuron.pos(:,1))*0.05; % cut off 10% on each side
edge12 = max(neuron.pos(:,1)) - edge11;
% edge21 = max(neuron.pos(:,2))*0.1; % cut off 10% on each side
% edge22 = max(neuron.pos(:,2)) - edge21;

binSize = 10;

% pos1 = 0:10:ceil(max(neuron.pos(:,1)));pos2 = 0:10:ceil(max(neuron.pos(:,2)));
% pos1 = edge11:binSize:edge12;pos2 = edge21:binSize:edge22;
pos1 = edge11:binSize:edge12;
idx1 = find(neuron.pos(:,1) >= edge11, 1, 'first');
idx2 = find(neuron.pos(:,1) <= edge12, 1, 'last');
pos2 = min(neuron.pos(idx1:idx2,2)):binSize:max(neuron.pos(idx1:idx2,2));


% [xpos1,ypos2] = meshgrid(pos1, pos2);
countTime = zeros(length(pos1),length(pos2));
ts = 1;
while ts < length(behav.time)
    %     ts
    [~,idxxi] = find(pos1 <= behav.position(ts,1), 1, 'last');
    [~,idyyi] = find(pos2 <= behav.position(ts,2), 1, 'last');
    for j = ts+1:length(behav.time)
        [~,idxxj] = find(pos1 <= behav.position(j,1), 1, 'last');
        [~,idyyj] = find(pos2 <= behav.position(j,2), 1, 'last');
        if idxxj == idxxi & idyyj == idyyi
            countTime(idxxi,idyyi) = countTime(idxxi,idyyi)+behav.time(j)-behav.time(j-1);
        else
            ts = j;
            break;
        end
    end
    if ts < j
        break;
    end
end
countTime = countTime'/1000;

% figure
% imagesc(countTime);
% % pcolor(xpos1,ypos2,countTime);
% % set(h,'EdgeColor','none')
% caxis ([0,max(countTime(:))])
% newmap = jet(256);
% newmap(1,:) = [1 1 1];
% colormap(newmap);
% colorbar
% figure
% countTime2 = nan(size(countTime)+1);
% countTime2(1:end-1,1:end-1) = countTime;
% countTime2(countTime2 == 0) = NaN;
% pcolor(countTime2);
% colormap(jet)
% colorbar;
% set(gca, 'color', 'w', 'ydir', 'reverse')
% shading flat;
% % shading interp
% axis image

%maxCount = zeros(1,num);
maxCount = 5; %add this line when you obtain the maximum of maxCount after running once
% firingRateAll = zeros(length(pos1),length(pos2),num);
% countAll = zeros(length(pos1),length(pos2),num);
firingRateAll = cell(1,num);
countAll = cell(1,num);
for k = 1:num
    
    if strcmpi(temp,'trace')
        if threshold < 1
            thresh = (max(neuron.trace(Segment(k),:))-min(neuron.trace(Segment(k),:)))*threshold; % the threshold above which the neuron is active
        else
            thresh = threshold(Segment(k));
        end
        idx = find(neuron.trace(Segment(k),:)>thresh);
    elseif strcmpi(temp,'S')
        if threshold < 1
            thresh = (max(neuron.trace(Segment(k),:))-min(neuron.trace(Segment(k),:)))*threshold; % the threshold above which the neuron is active
        else
            thresh = threshold(Segment(k));
        end
        idx = find(neuron.S(Segment(k),:)>thresh);
    end
    if ~isempty(idx)
        count = countingFiringBins(idx,neuron);
        firingRate = count./countTime;
        firingRate(countTime == 0) = NaN;
        
        firingRate2 = nan(size(firingRate)+1);
        firingRate2(1:end-1,1:end-1) = firingRate;
        countTime2 = nan(size(countTime)+1);
        countTime2(1:end-1,1:end-1) = countTime;
        firingRate2(countTime2 == 0) = NaN;
        %         maxCount = max(firingRate(:));
        if plotting
            figure;
            %              pcolor(firingRate2);
            firingRateSmoothing = filter2DMatrices(firingRate2, 1);
            pcolor(firingRateSmoothing);
            colormap(jet)
            %  maxCount(k) = max(firingRate(:)); % comment this line when you obtain the maximum of maxCount after running once
            caxis([0,max(maxCount)])
            colorbar;
            set(gca, 'color', 'w', 'ydir', 'reverse')
            shading flat;
            if intensity
                shading interp
            end
            axis image
            axis off equal;
            
            title(['Cell #', num2str(Segment(k))],'FontName','Arial','FontSize',10,'FontWeight','bold')
            hold off
            saveas(gcf,fullfile(fpath,strcat('CellSpatialMatch',num2str(Segment(k)),'.tif')))
            saveas(gcf,fullfile(fpath,strcat('CellSpatialMatch',num2str(Segment(k)),'.fig')))
        end
        
        firingRateAll{k} = firingRate;
        countAll{k} = count;
        
        %     firingRateAll(:,:,k) = firingRate;
        %     countAll(:,:,k) = count;
    end
    
    %     max(maxCount)
    
end

function count=countingFiringBins(idx,neuron)
global pos1 pos2
count = zeros(length(pos1),length(pos2));
for i = 1:length(idx)
    [~,idxx] = find(pos1 <= neuron.pos(idx(i),1), 1, 'last');
    [~,idyy] = find(pos2 <= neuron.pos(idx(i),2), 1, 'last');
    count(idxx,idyy) = count(idxx,idyy)+1;
end
count = count';
% count(end+1,:) = count(end,:);
% count(:,end+1) = count(:,end);

