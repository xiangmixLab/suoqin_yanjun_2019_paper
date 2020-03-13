function plottingFiringBehaviorSpatialForSingleData_adapted_dryland(neuron,behavpos,behavtime,behavROI,firingRate,segment,threshFiring,threshSpatial,threshSpatial_3,conditionfolder,sectionindex,binsize,tempname,object,firingRate3,labfr3,countTime,loc,j)

%% This function is used to overlay the neuron activity onto behaviors
experimentName=['start',num2str(loc)];

folderName = [conditionfolder{sectionindex},'/','FiguresFiringBehaviorSpatial'];
if ~exist(folderName,'dir')
    mkdir(folderName);
end
fpath=[folderName '/'];

for i=1:size(behavpos,2)
    downsampling = length(neuron{i}.time)/size(neuron{i}.C,2);
    if downsampling ~= 1
        neuron{i}.time = double(neuron{i}.time);
        neuron{i}.time = neuron{i}.time(1:downsampling:end);
    end
    temp = find(diff(behavtime{i})<=0);
    behavtime{i}(temp+1) = behavtime{i}(temp)+1;
    [ti,tiIdx]=unique(neuron{i}.time);
    neuron{i}.C = interp1(ti,neuron{i}.C(:,tiIdx)',behavtime{i})'; %%
    neuron{i}.S = interp1(ti,neuron{i}.S(:,tiIdx)',behavtime{i})'; %%
    neuron{i}.trace = interp1(ti,neuron{i}.trace(:,tiIdx)',behavtime{i})'; %%

end

plot_row=2;
kk=0;
numFig=10;
for i = segment
    disp(mod(kk,numFig))
    kk = kk+1;
    if mod(kk-1,numFig) == 0
        ax = figure('visible','off');
        set(ax, 'Position', [100, 100, 1000, 1000]);
        k = 0;
    end
    
    
    %   thresh = (max(neuron.S(i,:))-min(neuron.S(i,:)))*0.1; % the threshold above which the neuron is active
    thresh = threshFiring(i);
     k = k+1;
    for q=1:size(neuron,2)
        idx{q} = neuron{q}.S(i,:)>thresh;
        if isequal(tempname,'S')||isempty(tempname)
           idx{q} = neuron{q}.S(i,:)>thresh;
           max_pk{q}=max(neuron{q}.S(:));
        end
        if isequal(tempname,'trace')
           idx{q} = neuron{q}.trace(i,:)>thresh;
           max_pk{q}=max(neuron{q}.trace(:));
        end
        if isequal(tempname,'C')
           idx{q} = neuron{q}.C(i,:)>thresh;
           max_pk{q}=max(neuron{q}.C(:));
        end
    end
    
%% ploting raw trace and its firing
%     subplot(numFig,plot_row,plot_row*k-(plot_row-1))
%     plot(neuron.C(i,:), 'b')
%     hold on
%     title(['Cell' num2str(i)],'FontSize',8,'FontName','Arial')
%     set(gca,'FontSize',8)
%     axis tight
%     ylim([0 max_pk]);% temporally added
%     hold off
%     set(gca,'Xtick',[])
%     if mod(kk,numFig) == 0 || kk == max(segment)
%         set(gca,'Xtick',[1 ceil(neuron.num2read/2) neuron.num2read])
%     end
    
    %% plotting animal behavior trajectries    
    subplot(numFig,plot_row,plot_row*k-(plot_row-2))
    for q=1:size(neuron,2)
        plot(behavpos{q}(:,1),behavpos{q}(:,2),'k')
        hold on
        idx{q}=idx{q}(1:length(behavpos{q}(:,1)));
        plot(behavpos{q}(idx{q},1),behavpos{q}(idx{q},2),'r.','MarkerSize',1)
        plot(0,0);
        plot(behavROI(1,3),behavROI(1,4));
    end
    
    posObjects=object;
    if sum(posObjects)~=0
        for i5 = 1:size(posObjects,1)
            scatter(posObjects(i5,1),behavROI(1,4)-posObjects(i5,2)+1,binsize*2,'k','filled')
        end
    end

    title(['Cell' num2str(i)],'FontSize',8,'FontName','Arial')
    set(gca,'FontSize',8)
    axis image
     axis ij
%         xlim([min(neuron.pos(:,1)) max(neuron.pos(:,1))]);
%         ylim([min(neuron.pos(:,2)) max(neuron.pos(:,2))]);
    %     plot(ms.pos(idx2,1),ms.pos(idx2,2),'r.')
    hold off
    
     %% plot cell firing map
%     if isempty(firingRate{i})
%         for iii=1:length(firingRate)
%             if ~isempty(firingRate{1,iii})          
%                 [mm,nn]=size(firingRate{1,iii});
%                 break;
%             end
%         end
%         firingRate{i}=zeros(mm,nn);
%     end
%     ft1=firingRate{i};
%     ft1(isnan(ft1))=0;
%     ft1(ft1==inf)=0;
% 
%     firingRateSmoothing = filter2DMatrices(ft1, 1);
%     firingRateSmoothing(countTime==0)=nan;
%     firingRateSmoothing2 = nan(size(firingRateSmoothing)+1);
%     firingRateSmoothing2(1:end-1,1:end-1) = firingRateSmoothing;
%     hand=subplot(numFig,plot_row,plot_row*k-(plot_row-3));
% 
%     try
%         pcolor(firingRateSmoothing2);
%         hold on;
%         colormap(jet)
%         caxis([0,max(firingRateSmoothing2(:))]);
%         axis ij;
%         axis image
%         posObjects=round(object./binsize);
%         if sum(posObjects)~=0
%             for i5 = 1:size(posObjects,1)
%                 scatter(posObjects(i5,1)+1,size(firingRateSmoothing,1)-posObjects(i5,2)+1,binsize*2,'k','filled')
%             end
%         end
%         
%         shading flat;
%         title(['event count'],'FontName','Arial','FontSize',8,'FontWeight','bold')
%         hold off
%         if mod(kk,numFig) ~= 0 
%             set(gca,'Xtick',[])
%         end        
%         colorbar('position',[0.65,hand.Position(2),0.01,hand.Position(4)]);
%     catch
%         continue;
%     end
% 
%     %% plot cell firing map 3 (firing rate by default, change the parameter outside to change this)
%     if isempty(firingRate3{i})
%         for iii=1:length(firingRate3)
%             if ~isempty(firingRate3{1,iii})          
%                 [mm,nn]=size(firingRate3{1,iii});
%                 break;
%             end
%         end
%         firingRate3{i}=zeros(mm,nn);
%     end
%     ft3=firingRate3{i};
%     ft3(isnan(ft3))=nan;
%     ft3(ft3==inf)=nan;
% 
%     firingRateSmoothing = filter2DMatrices(ft3, 1);
%     firingRateSmoothing(countTime==0)=nan;
%     firingRateSmoothing2 = nan(size(firingRateSmoothing)+1);
%     firingRateSmoothing2(1:end-1,1:end-1) = firingRateSmoothing;
%     hand=subplot(numFig,plot_row,plot_row*k-(plot_row-4));
% 
%     try
%         pcolor(firingRateSmoothing2);
%         hold on;
%         colormap(jet)
%         caxis([0,max(firingRateSmoothing2(:))]);
%         axis ij;
%         axis image
%         posObjects=round(object./binsize);
%         if sum(posObjects)~=0
%             for i5 = 1:size(posObjects,1)
%                 scatter(posObjects(i5,1)+1,size(firingRateSmoothing,1)-posObjects(i5,2)+1,binsize*2,'k','filled')
%             end
%         end
%         
%         shading flat;
%         hold off
% 
%         colorbar('position',[0.9,hand.Position(2),0.01,hand.Position(4)]);
%         
%         if mod(kk,numFig) ~= 0 
%             set(gca,'Xtick',[])
%         end
% 
%     catch
%         continue;
%     end
    set(gcf,'renderer','painters');
    if mod(kk,numFig) == 0 || kk == max(segment)
        saveas(gcf,fullfile(fpath,strcat(experimentName,'_CellFiringBeaviorSpatial_',num2str(i),'_binsize',num2str(binsize),tempname,'_','obj',num2str(j),'.fig')))
        saveas(gcf,fullfile(fpath,strcat(experimentName,'_CellFiringBeaviorSpatial_',num2str(i),'_binsize',num2str(binsize),tempname,'_','obj',num2str(j),'.tif')))
        saveas(gcf,fullfile(fpath,strcat(experimentName,'_CellFiringBeaviorSpatial_',num2str(i),'_binsize',num2str(binsize),tempname,'_','obj',num2str(j),'.eps')),'epsc')
        close all;
        disp(['fin ',strcat(experimentName,'_CellFiringBeaviorSpatial_',num2str(i),'_binsize',num2str(binsize),tempname,'_','obj',num2str(j),'.fig')])
    end   
end
