function [firingRateAll,countAll,countTimeAll,countTime,amplitudeAll,binInfo] = calculatingCellSpatialForSingleData(neuron,behav,binsize,segments,threshold,temp,intensity,plotting,countTimeThresh)
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
if ~exist('plotting','var') || isempty(plotting)
    plotting = false;
end
if ~exist('intensity','var') || isempty(intensity)
    intensity = false;
end
if ~exist('temp','var') || isempty(temp)
    temp = 'S';
end
if ~exist('threshold','var') || isempty(threshold)
    threshold = 0.1;
end
if ~exist('segments','var') || isempty(segments)
    segments = 1:size(neuron.C,1);
end
if ~exist('binsize','var') || isempty(binsize)
    binsize = 15;%attention: value has been 5, 10, 15
end
if ~exist('countTimeThresh','var') || isempty(countTimeThresh)
    countTimeThresh = [0.2 10]; % unit: sec
end

behavpos = behav.position;behavtime = behav.time;behavROI = behav.ROI;

downsampling = length(neuron.time)/size(neuron.trace,2);
if downsampling ~= 1
    %     downsampling == 2
    neuron.time = double(neuron.time);
    neuron.time = neuron.time(1:downsampling:end);
end
t = find(diff(behavtime)<=0);
while ~isempty(t)
    behavtime(t+1) = behavtime(t)+1;
    t = find(diff(behavtime)<=0);
end
if isempty(neuron.pos)
neuron.pos = interp1(behavtime,behavpos,neuron.time); %%
end


folderName = 'FiguresCellSpatial';
if ~exist(folderName,'dir')
    mkdir(folderName)
end

fpath=[folderName];
global ts pos1 pos2
num = length(segments);

pos1 = 0:binsize:ceil(max(neuron.pos(:,1)));
pos2 = 0:binsize:ceil(max(neuron.pos(:,2)));
% pos1 = 0:binsize:round(max(neuron.pos(:,1)));
% pos2 = 0:binsize:round(max(neuron.pos(:,2)));
% pos1 = 0:binsize:round(max(neuron.pos(:,1)));
% pos2 = 0:binsize:round(max(neuron.pos(:,2)));
% pos1 = floor(min(neuron.pos(:,1))):binsize:ceil(max(neuron.pos(:,1)));
% pos2 = floor(min(neuron.pos(:,2))):binsize:ceil(max(neuron.pos(:,2)));
% decrement1 = max(neuron.pos(:,1)-min(neuron.pos(:,1)))*0.005;
% decrement2 = max(neuron.pos(:,2)-min(neuron.pos(:,2)))*0.005;
% pos1 = floor(min(neuron.pos(:,1))+decrement1):binsize:ceil(max(neuron.pos(:,1))-decrement1);
% pos2 = floor(min(neuron.pos(:,2))+decrement2):binsize:ceil(max(neuron.pos(:,2))-decrement2);

% sprintf("The number of bins is %dx%d",length(pos1),length(pos2))

% if max(pos1) < floor(max(neuron.pos(:,1)))
%     pos1 = [pos1 max(pos1)+binsize];
% end
% if max(pos2) < floor(max(neuron.pos(:,2)))
%     pos2 = [pos2 max(pos2)+binsize];
% end
binInfo.binsize = binsize;binInfo.pos1 = pos1;binInfo.pos2 = pos2;
binInfo.xpos = [min(neuron.pos(:,1)),max(neuron.pos(:,1))];
binInfo.ypos = [min(neuron.pos(:,2)),max(neuron.pos(:,2))];
% pos1 = 0:binsize:ceil(behavROI(:,3));%we changed to ROI
% pos2 = 0:binsize:ceil(behavROI(:,4));

countTime = zeros(length(pos1),length(pos2));
ts = 1;
while ts < length(behavtime)
    %     ts
    [~,idxxi] = find(pos1 <= behavpos(ts,1), 1, 'last');
    [~,idyyi] = find(pos2 <= behavpos(ts,2), 1, 'last');
    for j = ts+1:length(behavtime)
        [~,idxxj] = find(pos1 <= behavpos(j,1), 1, 'last');
        [~,idyyj] = find(pos2 <= behavpos(j,2), 1, 'last');
        if idxxj == idxxi & idyyj == idyyi
            countTime(idxxi,idyyi) = countTime(idxxi,idyyi)+behavtime(j)-behavtime(j-1);
        else
            ts = j;
            break;
        end
    end
    if ts < j
        break;
    end
end
countTime = countTime'/1000;%purpose: because behavtime is recorded in milisec.
% countTime = countTime';
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
countTimeAll = cell(1,num);
amplitudeAll = cell(1,num);
for k = 1:num
    if strcmpi(temp,'trace')
        if length(threshold) <= 1
            %         thresh = (max(neuron.trace(Segment(k),:))-min(neuron.trace(Segment(k),:)))*threshold; % the threshold above which the neuron is active
            thresh = (max(neuron.trace(segments(k),:))-0)*threshold; % the threshold above which the neuron is active
        else
            thresh = threshold(segments(k));
        end
        %         idx = find(neuron.trace(segments(k),:)>thresh);
        [pks,locs] = findpeaks(neuron.trace(segments(k),:),'MinPeakHeight',thresh);
        idx = locs;
    elseif strcmpi(temp,'S')
        if length(threshold) <= 1
            %         thresh = (max(neuron.S(Segment(k),:))-min(neuron.S(Segment(k),:)))*threshold; % the threshold above which the neuron is active
            thresh = (max(neuron.S(segments(k),:))-0)*threshold; % the threshold above which the neuron is active
        else
            thresh = threshold(segments(k));
        end
        idx = find(neuron.S(segments(k),:)>thresh);
    elseif strcmpi(temp,'C')
        if length(threshold) <= 1
            %         thresh = (max(neuron.C(Segment(k),:))-min(neuron.C(Segment(k),:)))*threshold; % the threshold above which the neuron is active
            thresh = (max(neuron.C(segments(k),:))-0)*threshold; % the threshold above which the neuron is active
        else
            thresh = threshold(segments(k));
        end
        %         idx = find(neuron.C(segments(k),:)>thresh);
        %         [pks,locs,w,p] = findpeaks(neuron.C(segments(k),:),'MinPeakDistance',100,'MinPeakHeight',thresh);
        [pks,locs] = findpeaks(neuron.C(segments(k),:),'MinPeakHeight',thresh);
        idx = locs;
    end
    if ~isempty(idx)
        count = countingFiringBins(idx,neuron);
      %  amplitude = countingAmplitudeFiringBins(idx,neuron,segments(k));
        %% 0.2sec threshold && 10sec threshold
%         countTime_smaller_than_thr=countTime<countTimeThresh(1); %0.2sec
%         count(countTime_smaller_than_thr)=0;
%         countTime(countTime_smaller_than_thr)=0;
      %  amplitude = amplitude./count;amplitude(countTime == 0) = NaN;
      amplitude = [];
        
        %         countTime_larger_than_thr=countTime>countTimeThresh(2); %10sec
        %         count(countTime_larger_than_thr)=0;
        %         countTime(countTime_larger_than_thr)=0;
        
        firingRate = count./countTime;
        firingRate(countTime == 0) = NaN;
        %  firingRate(isnan(firingRate))=0;
        
        firingRate2 = nan(size(firingRate)+1);
        firingRate2(1:end-1,1:end-1) = firingRate;
        countTime2 = nan(size(countTime)+1);
        countTime2(1:end-1,1:end-1) = countTime;
        firingRate2(countTime2 == 0) = NaN;
        %  firingRate2(isnan(firingRate2))=0;
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
            axis off;
            axis ij
            title(['Cell #', num2str(segments(k))],'FontName','Arial','FontSize',10,'FontWeight','bold')
            hold off
            saveas(gcf,fullfile(fpath,['CellSpatialMatch',num2str(segments(k)),'.tif']))
            saveas(gcf,fullfile(fpath,['CellSpatialMatch',num2str(segments(k)),'.fig']))
        end
        
        firingRateAll{k} = firingRate;
        countAll{k} = count;
        countTimeAll{k} = countTime;
        amplitudeAll{k} = amplitude;
        %     firingRateAll(:,:,k) = firingRate;
        %     countAll(:,:,k) = count;
    end
    
    %max(maxCount)
    
end

function count=countingFiringBins(idx,neuron)
global pos1 pos2
count = zeros(length(pos1),length(pos2));
for i = 1:length(idx)
    % if idx(i)<=length(neuron.pos) %sometimes in behav data, the position is less longer than neuron.S, which cause crash
    [~,idxx] = find(pos1 <= neuron.pos(idx(i),1), 1, 'last');
    [~,idyy] = find(pos2 <= neuron.pos(idx(i),2), 1, 'last');
    count(idxx,idyy) = count(idxx,idyy)+1;
    % end
end
count = count';
% count(end+1,:) = count(end,:);
% count(:,end+1) = count(:,end);

function amplitude = countingAmplitudeFiringBins(idx,neuron,k)
global pos1 pos2
amplitude = zeros(length(pos1),length(pos2));
for i = 1:length(idx)
    % if idx(i)<=length(neuron.pos) %sometimes in behav data, the position is less longer than neuron.S, which cause crash
    [~,idxx] = find(pos1 <= neuron.pos(idx(i),1), 1, 'last');
    [~,idyy] = find(pos2 <= neuron.pos(idx(i),2), 1, 'last');
    amplitude(idxx,idyy) = amplitude(idxx,idyy)+neuron.S(k,idx(i));
    % end
end
amplitude = amplitude';

