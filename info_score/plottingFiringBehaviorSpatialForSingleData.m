function plottingFiringBehaviorSpatialForSingleData(neuron,behav,firingRate,segments,threshFiring,temp,threshSpatial,experimentName)

%% This function is used to overlay the neuron activity onto behaviors
% Inputs:
%        (1) neuron: a source2D variable, including identified neurons with traces and spatial information, which is obtained by runing cnmfe codes
%        (2) behav: behavior information, which is obtained by using Tristan's code
%        (3) Segment: a vector, e.g,1:10 (display the traces of the first 10 identified neurons)
%        (4) downsampling rate [e.g., by a factor of 2]

% Important parameter in the code: thresh, the threshold above which the neuron is active. By default, it is 10% of the maximum trace value of each neuron
if ~exist('temp','var') || isempty(temp)
    temp = 'S';
end
if strcmpi(temp,'trace')
    dataFiring = neuron.trace;
elseif strcmpi(temp,'S')
    dataFiring = neuron.S;
elseif strcmpi(temp,'C')
    dataFiring = neuron.C;
end


folderName = fullfile('results','figures','FiringBehaviorSpatial');
if ~exist(folderName, 'dir')
    mkdir(folderName);
end

fpath=folderName;

% figure

downsampling = length(neuron.time)/size(neuron.trace,2);
if downsampling ~= 1
    %     downsampling == 2
    neuron.time = double(neuron.time);
    neuron.time = neuron.time(1:downsampling:end);
end
temp0 = find(diff(behav.time)<=0);
behav.time(temp0+1) = behav.time(temp0)+1;
neuron.pos = interp1(behav.time,behav.position,neuron.time); %%

%
% thresh = 0;
numFig = 10;
k = 0;kk = 0;
for t = 1:length(segments)
    i = segments(t);
    kk = kk+1;
    if mod(kk-1,numFig) == 0
        ax = figure;
        set(ax, 'Position', [100, 100, 600, 1000]);
        k = 0;
        ha = tight_subplot(numFig,3,[.015 .02],[.03 .02],[.1 .1]);
    end
    
    
    %   thresh = (max(neuron.S(i,:))-min(neuron.S(i,:)))*0.1; % the threshold above which the neuron is active
    thresh = threshFiring(i);
    k = k+1;
    if strcmpi(temp,'S')
        idx = dataFiring(i,:)>thresh;
    elseif strcmpi(temp,'C') | strcmpi(temp,'trace')
        [~,locs] = findpeaks(dataFiring(i,:),'MinPeakHeight',thresh);
        idx = locs;
    end
    %% ploting raw trace and its firing
    %     subplot(length(segment),3,3*k-2)
    %     subplot(numFig,3,3*k-2)
    axes(ha(3*k-2))
    %     subplot(5,2,plotPositionFiring)
    plot(neuron.C(i,:), 'b')
    hold on
    plot(neuron.S(i,:),'r')
    line(get(gca,'XLim'),[thresh thresh],'Color','k','LineStyle','--');
    title(['Cell #', num2str(i)],'FontSize',8,'FontName','Arial')
    set(gca,'FontSize',8)
    %     plot([0 neuron.num2read],thresh*[1 1],'k--')
    axis tight
    %     xlim([0 neuron.num2read]);
    %     ylim([min([neuron.trace(i,:);neuron.firing(:,i)]) max([neuron.trace(:,i);neuron.firing(:,i)])])
    % ylim([min([ms.trace(:,i)]) max([ms.trace(:,i)])])
    hold off
    set(gca,'Xtick',[])
    if mod(kk,numFig) == 0 || kk == max(segments)
        set(gca,'Xtick',[1 ceil(neuron.num2read/2) neuron.num2read])
    end
    
    %% plotting animal behavior trajectries
    %     subplot(length(segment),3,3*k-1)
    %  subplot(numFig,3,3*k-1)
    axes(ha(3*k-1))
    %     subplot(5,2,plotPositionPos)
    plot(neuron.pos(:,1),neuron.pos(:,2),'k')
    hold on
    plot(neuron.pos(idx,1),neuron.pos(idx,2),'r.')
    
    if isfield(behav,'object') & sum(behav.object(:)) ~=0
        scatter((behav.object(:,1)-behav.ROI(1))*behav.trackLength/behav.ROI(3),max(behav.position(:,2))-behav.object(:,2)*behav.trackLength/behav.ROI(3),60,'k','filled')
    end
    
    title(['Cell #' num2str(i)],'FontSize',8,'FontName','Arial')
    set(gca,'FontSize',8)
    axis image
    axis ij
    %     xlim([min(neuron.pos(:,1)) max(neuron.pos(:,1))]);
    %     ylim([min(neuron.pos(:,2)) max(neuron.pos(:,2))]);
    %     plot(ms.pos(idx2,1),ms.pos(idx2,2),'r.')
    hold off
    set(gca,'Xtick',[])
    
    %% plot cell firing map
    maxCount = threshSpatial;
    firingRate_trim = firingRate{i}(:,~all(isnan(firingRate{i})));
    firingRate_trim = firingRate_trim(~all(isnan(firingRate_trim),2),:);
    firingRate_trim(isnan(firingRate_trim)) = 0;
    firingRateSmoothing = filter2DMatrices(firingRate_trim, 1);
    
%     
%     firingRateSmoothing = filter2DMatrices(firingRate{i}, 1);
    firingRateSmoothing2 = nan(size(firingRateSmoothing)+1);
    firingRateSmoothing2(1:end-1,1:end-1) = firingRateSmoothing;
    %      subplot(length(segment),3,3*k)
    % subplot(numFig,3,3*k)
    axes(ha(3*k))
    try
        pcolor(firingRateSmoothing2);
        colormap(jet)
        caxis([0,max(maxCount)])
        %         colorbar('eastoutside');
        % set(gca, 'color', 'w', 'ydir', 'reverse')
        shading flat;
        %         axis square
        axis image
        axis ij
        % axis off
        %                 if sum(behav.object(:)) ~=0
        %          scatter((behav.object(:,1)-behav.ROI(1))*behav.trackLength/behav.ROI(3),max(behav.position(:,2))-behav.object(:,2)*behav.trackLength/behav.ROI(3),60,'k','filled')
        %                 end
        
        title(['Cell ', num2str(i),' (peak FR: ',num2str(max(firingRate{i}(:))),')'],'FontName','Arial','FontSize',8,'FontWeight','bold')
        hold off
        
        if mod(kk,numFig) ~= 0
            set(gca,'Xtick',[])
        end
        if mod(kk,numFig) == 0
            %    colorbar('eastoutside');
            c = colorbar;
            % c.Position = [0.81 .4 .015 .2];
            c.Position = [0.9 0.4 .015 .1];
            % c.Location = 'east';
            c.Label.String = 'smoothed FR (Hz)';
            c.Label.FontSize = 8;%c.Label.FontWeight = 'bold';
            c.FontSize = 8;
            
        end
    catch
        continue;
    end
    
    if mod(kk,numFig) == 0 || kk == max(segments)
        saveas(gcf,fullfile(fpath,strcat(experimentName,'CellFiringBeaviorSpatial',num2str(i),'.fig')))
        saveas(gcf,fullfile(fpath,strcat(experimentName,'CellFiringBeaviorSpatial',num2str(i),'.png')))
    end
    
end
%     saveas(gcf,fullfile(fpath,strcat('CellFiringBeaviorSpatial','.tif')))
%     saveas(gcf,fullfile(fpath,strcat('CellFiringBeaviorSpatial','.fig')))

