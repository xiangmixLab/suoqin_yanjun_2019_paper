function dryland_single_mouse_illustration_forgeR(behavname,startlocation,savedir,fileloc)

for i=1:length(startlocation)
    savedir1{i}=[savedir,'\','trial',num2str(i),'_','beforeObj'];
end
for i=1:length(startlocation)
    savedir1{i+length(startlocation)}=[savedir,'\','trial',num2str(i),'_','AfterObj'];
end
for i=1:length(startlocation)
    savedir1{i+2*length(startlocation)}=[savedir,'\','trial',num2str(i),'_','allObj'];
end

for i=1:length(savedir1)
    mkdir(savedir1{i});
end
%% behav track split
originalPos={};
originalTime={};
posBeforeObj={};
timeBeforeObj={};
posAfterObj={};
timeAfterObj={};
dis_to_obj_all={};
reach_obj_point_all={};
away_from_obj_point_all={};
ROIlist={};
objects={};
objmodifier=[0 -0];%amend position

for i=1:length(behavname)
    load(behavname{i});
    object=behav.object;
    object(:,1)=object(:,1)+objmodifier(1);
    object(:,2)=behav.ROI(4)-(object(:,2)+objmodifier(2));
    behavpos=behav.position;
    behavtime=behav.time;
    
    dis_to_obj={};
    reach_obj_point={};
    away_from_obj_point={};
    for j=1:size(object,1)
        coor_diff=(behavpos-repmat(object(j,:),size(behavpos,1),1));
        dis_to_obj{j}=(coor_diff(:,1).^2+coor_diff(:,2).^2).^0.5;
        dis_to_obj{j}(isnan(dis_to_obj{j}))=inf;
        reach_obj_point{j}=min(find(dis_to_obj{j}<10));
        if isempty(reach_obj_point{j})
            reach_obj_point{j}=min(find(dis_to_obj{j}==min(dis_to_obj{j})));
        end
        if ~isempty(min(find(dis_to_obj{j}(reach_obj_point{j}:end)>40)))
            away_from_obj_point{j}=min(find(dis_to_obj{j}(reach_obj_point{j}:end)>40))+reach_obj_point{j};
        else
            away_from_obj_point{j}=length(behavtime);
        end
        posBeforeObj{i,j}=behavpos(1:reach_obj_point{j},:);
        timeBeforeObj{i,j}=behavtime(1:reach_obj_point{j});
        posAfterObj{i,j}=behavpos(reach_obj_point{j}:away_from_obj_point{j},:);
        timeAfterObj{i,j}=behavtime(reach_obj_point{j}:away_from_obj_point{j});
    end

    originalPos{i}=behavpos;
    originalTime{i}=behavtime;  
    dis_to_obj_all{i}=dis_to_obj;
    reach_obj_point_all{i}=reach_obj_point;
    away_from_obj_point_all{i}=away_from_obj_point;
    ROIlist{i}=behav.ROI;
    objects{i}=object;
end

save([savedir,'\','behavTrackData.mat'],'originalPos','originalTime','posBeforeObj','timeBeforeObj','posAfterObj','timeAfterObj','dis_to_obj_all','reach_obj_point_all','away_from_obj_point_all','ROIlist','objects');

%% background frame
    load(behavname{4});
    count=1;
    frame1=[];
    try
        for i=1:10:behav.numFrames
            frame1(:,:,:,count)=double(msReadFrame(behav,i,false,false,false))/255;
            count=count+1;
        end
        frame_bg=median(frame1,4);
    catch
        frame_bg=imread([savedir,'\bg.tif']);
    end

%% combine behav track plot, before reach obj
    
    for loc=unique(startlocation)
        figure;
        imshow(frame_bg);
        hold on;
        for i=1:length(startlocation)
            if startlocation(i)==loc
                colorTraj = distinguishable_colors(size(posBeforeObj,2));
                for j=1:size(objects{i},1)
                    plot((posBeforeObj{i,j}(:,1)+behav.ROI(1))*behav.ROI3/behav.trackLength,(posBeforeObj{i,j}(:,2)+behav.ROI(2))*behav.ROI3/behav.trackLength,'-','color',colorTraj(i,:));
    %                 drawnow;
    %                 pause(0.5);
                    firstLoc=min(find(~isnan(posBeforeObj{i,j}(:,1))));
                    plot((posBeforeObj{i,j}(firstLoc,1)+behav.ROI(1))*behav.ROI3/behav.trackLength,(posBeforeObj{i,j}(firstLoc,2)+behav.ROI(2))*behav.ROI3/behav.trackLength,'c.','markerSize',20);
%                     text((posBeforeObj{i,j}(firstLoc,1)+behav.ROI(1))*behav.ROI3/behav.trackLength+2,(posBeforeObj{i,j}(firstLoc,2)+behav.ROI(2))*behav.ROI3/behav.trackLength+2,num2str(i),'color','w');
     %                 plot((originalPos{i}(:,1)+behav.ROI(1))*behav.ROI3/behav.trackLength,(originalPos{i}(:,2)+behav.ROI(2))*behav.ROI3/behav.trackLength,'-','color',colorTraj(i,:));            
                    plot((objects{i}(j,1)+behav.ROI(1))*behav.ROI3/behav.trackLength,(objects{i}(j,2)+behav.ROI(2))*behav.ROI3/behav.trackLength,'c.','markerSize',40)
                end
            end
        end
        saveas(gcf,[savedir1{find(unique(startlocation)==loc)},'\','traj_before_reach_obj_start_',num2str(loc),'.fig']);
        saveas(gcf,[savedir1{find(unique(startlocation)==loc)},'\','traj_before_reach_obj_start_',num2str(loc),'.eps'],'epsc');
        saveas(gcf,[savedir1{find(unique(startlocation)==loc)},'\','traj_before_reach_obj_start_',num2str(loc),'.tif']);
    end
    
%% combine behav track plot and neuron, before reach obj
    load([fileloc,'\','neuronIndividuals_new.mat']);
    load([fileloc,'\','further_processed_neuron_extraction_final_result.mat']);
    
    temp='S';
    
%     maxS = max(neuron.C,[],2);
    thresh =3*std(neuron.S,[],2);
    binsize=5;

    for loc=unique(startlocation)
        firingrateAll={};
        countAll={};
        countTimeAll={};
        amplitudeAll={};
        countt=1;
        for i=1:length(startlocation)
            if startlocation(i)==loc
                object1=objects{i};
                object1(:,2)=ROIlist{i}(4)-(object1(:,2));

                for j=1:size(object1,1)
                    n=Sources2D();
                    n.C=neuronIndividuals_new{i}.C(:,1:reach_obj_point_all{i}{j}/2);
                    n.S=neuronIndividuals_new{i}.S(:,1:reach_obj_point_all{i}{j}/2);
                    n.trace=neuronIndividuals_new{i}.trace(:,1:reach_obj_point_all{i}{j}/2);
                    n.time=neuronIndividuals_new{i}.time(1:reach_obj_point_all{i}{j});

                    [firingrate,count,~,countTime] = calculatingCellSpatialForSingleData_Suoqin(n,posBeforeObj{i,j},timeBeforeObj{i,j},ROIlist{i},binsize,1:size(n.C,1),thresh,temp,[],[],[0 1000000]);
                    [~,~,~,~,amplitude,amplitude_normalized] = calculatingCellSpatialForSingleData_Suoqin(n,posBeforeObj{i,j},timeBeforeObj{i,j},ROIlist{i},binsize,1:size(n.C,1),thresh,temp,[],[],[0 1000000]);  %%%bin size suggests to be 15
    %                 [firingrate,count,~,countTime] = calculatingCellSpatialForSingleData_Suoqin(neuronIndividuals_new{i},originalPos{i},originalTime{i},ROIlist{i},binsize,1:size(neuronIndividuals_new{i}.C,1),thresh,temp,[],[],[0 1000000]);
    %                 [~,~,~,~,amplitude,amplitude_normalized] = calculatingCellSpatialForSingleData_Suoqin(neuronIndividuals_new{i},originalPos{i},originalTime{i},ROIlist{i},binsize,1:size(neuronIndividuals_new{i}.C,1),thresh,temp,[],[],[0 1000000]);  %%%bin size suggests to be 15
    %                 
                    for ipp=0:length(firingrate)-1
                        firingrateAll{countt+ipp,j}=firingrate{ipp+1};
                        countAll{countt+ipp,j}=count{ipp+1};
                        countTimeAll{countt+ipp,j}=countTime;
                        amplitudeAll{countt+ipp,j}=amplitude{ipp+1};
                        amplitudeNAll{countt+ipp,j}=amplitude_normalized{ipp+1};
                    end                                         
                end
                countt=countt+length(firingrate);
            end
        end
        
        for i=1:size(countTimeAll,1)
            size1(i)=size(countTimeAll{i,1},1);
            size2(i)=size(countTimeAll{i,1},2);
        end
        size1m=max(size1);
        size2m=max(size2);
        countTime={};
        for j=1:size(countTimeAll,2)
            countTime{j}=zeros(size1m,size2m);
        end
        for i=1:size(countTimeAll,1)
            for j=1:size(countTimeAll,2)
                if ~isempty(countTimeAll{i,j})
                    countTime{j}=countTime{j}+imresize(countTimeAll{i,j},[size1m,size2m]);
                end
                if ~isempty(firingrateAll{i,j})
                    firingrateAll{i,j}=imresize(firingrateAll{i,j},[size1m,size2m]);
                end
                if ~isempty(countAll{i,j})
                    countAll{i,j}=imresize(countAll{i,j},[size1m,size2m]);
                end
                if ~isempty(amplitudeAll{i,j})
                    amplitudeAll{i,j}=imresize(amplitudeAll{i,j},[size1m,size2m]);
                end
                if ~isempty(amplitudeNAll{i,j})
                    amplitudeNAll{i,j}=imresize(amplitudeNAll{i,j},[size1m,size2m]);
                end                
            end
        end        
        
        for j=1:size(countTimeAll,2)
%             [firingRateSmoothing1,sumFiringRateObject1,firingRateSmoothingt,~,radius,posObjects1,sumFiringRateObject_individual,~,~,sumFiringRateObject_2nd]=comparingFiringRateSingleConditionMultiObjects(countAll(:,j),binsize,countTime{j}, object1,'events',[],1,1,{'obj'});
%             save([savedir1{find(unique(startlocation)==loc)},'\','neuron_comparingCountevents_binsize',num2str(binsize),'_',temp,'trial',num2str(loc),'_obj',num2str(j),'.mat'],'firingRateSmoothing1','sumFiringRateObject1','posObjects1','countAll','object1','sumFiringRateObject_2nd');
%             saveas(gcf,[savedir1{find(unique(startlocation)==loc)},'\','neuron comparingCountevents_binsize',num2str(binsize),'_',temp,'trial',num2str(loc),'_obj',num2str(j),'.fig'],'fig');
%             saveas(gcf,[savedir1{find(unique(startlocation)==loc)},'\','neuron comparingCountevents_binsize',num2str(binsize),'_',temp,'trial',num2str(loc),'_obj',num2str(j),'.tif'],'tif');
%             saveas(gcf,[savedir1{find(unique(startlocation)==loc)},'\','neuron comparingCountevents_binsize',num2str(binsize),'_',temp,'trial',num2str(loc),'_obj',num2str(j),'.eps'],'epsc');
            [firingRateSmoothing2,sumFiringRateObject2,firingRateSmoothingt,~,radius,posObjects2,sumFiringRateObject_individual,~,~,sumFiringRateObject_2nd]=comparingFiringRateSingleConditionMultiObjects(countTimeAll(:,j)',binsize,countTime{j}, object1,'count Time',[],1,1,{'1','2','3','4','5','6','7','8','9','10'});
            save([savedir1{find(unique(startlocation)==loc)},'\','neuron_comparingCountTime_binsize',num2str(binsize),'_',temp,'trial',num2str(loc),'_obj',num2str(j),'_data.mat'],'firingRateSmoothing2','sumFiringRateObject2','posObjects2','countTimeAll','object1','sumFiringRateObject_2nd');
            saveas(gcf,[savedir1{find(unique(startlocation)==loc)},'\','neuron comparingCountTime_binsize',num2str(binsize),'_',temp,'trial',num2str(loc),'_obj',num2str(j),'.fig'],'fig');
            saveas(gcf,[savedir1{find(unique(startlocation)==loc)},'\','neuron comparingCountTime_binsize',num2str(binsize),'_',temp,'trial',num2str(loc),'_obj',num2str(j),'.tif'],'tif');
            saveas(gcf,[savedir1{find(unique(startlocation)==loc)},'\','neuron comparingCountTime_binsize',num2str(binsize),'_',temp,'trial',num2str(loc),'_obj',num2str(j),'.eps'],'epsc');
            [firingRateSmoothing3,sumFiringRateObject3,firingRateSmoothingt,~,radius,posObjects3,sumFiringRateObject_individual,~,~,sumFiringRateObject_2nd]=comparingFiringRateSingleConditionMultiObjects(firingrateAll(:,j)',binsize,countTime{j}, object1,'firing rate',[],1,1,{'1','2','3','4','5','6','7','8','9','10'});
            save([savedir1{find(unique(startlocation)==loc)},'\','neuron_comparingFiringRate_binsize',num2str(binsize),'_',temp,'trial',num2str(loc),'_obj',num2str(j),'_data.mat'],'firingRateSmoothing3','sumFiringRateObject3','posObjects3','firingrateAll','object1','sumFiringRateObject_2nd');
            saveas(gcf,[savedir1{find(unique(startlocation)==loc)},'\','neuron comparingFiringRate_binsize',num2str(binsize),'_',temp,'trial',num2str(loc),'_obj',num2str(j),'.fig'],'fig');
            saveas(gcf,[savedir1{find(unique(startlocation)==loc)},'\','neuron comparingFiringRate_binsize',num2str(binsize),'_',temp,'trial',num2str(loc),'_obj',num2str(j),'.tif'],'tif');
            saveas(gcf,[savedir1{find(unique(startlocation)==loc)},'\','neuron comparingFiringRate_binsize',num2str(binsize),'_',temp,'trial',num2str(loc),'_obj',num2str(j),'.eps'],'epsc');      
%             [firingRateSmoothing4,sumFiringRateObject4,firingRateSmoothingt,~,radius,posObjects4,sumFiringRateObject_individual,~,~,sumFiringRateObject_2nd]=comparingFiringRateSingleConditionMultiObjects(amplitudeAll(:,j),binsize,countTime{j}, object1,'Amplitude',[],1,1,{'obj'});
%             save([savedir1{find(unique(startlocation)==loc)},'\','neuron_comparingAmplitude_binsize',num2str(binsize),'_',temp,'trial',num2str(loc),'_obj',num2str(j),'_data.mat'],'firingRateSmoothing4','sumFiringRateObject4','posObjects4','amplitudeAll','object1','sumFiringRateObject_2nd');
%             saveas(gcf,[savedir1{find(unique(startlocation)==loc)},'\','neuron comparingAmplitude_binsize',num2str(binsize),'_',temp,'trial',num2str(loc),'_obj',num2str(j),'.fig'],'fig');
%             saveas(gcf,[savedir1{find(unique(startlocation)==loc)},'\','neuron comparingAmplitude_binsize',num2str(binsize),'_',temp,'trial',num2str(loc),'_obj',num2str(j),'.tif'],'tif');
%             saveas(gcf,[savedir1{find(unique(startlocation)==loc)},'\','neuron comparingAmplitude_binsize',num2str(binsize),'_',temp,'trial',num2str(loc),'_obj',num2str(j),'.eps'],'epsc');
            [firingRateSmoothing5,sumFiringRateObject5,firingRateSmoothingt,~,radius,posObjects5,sumFiringRateObject_individual,~,~,sumFiringRateObject_2nd]=comparingFiringRateSingleConditionMultiObjects(amplitudeNAll(:,j)',binsize,countTime{j}, object1,'Amplitude',[],1,1,{'1','2','3','4','5','6','7','8','9','10'});
            save([savedir1{find(unique(startlocation)==loc)},'\','neuron_comparingAmplitudeNorm_binsize',num2str(binsize),'_',temp,'trial',num2str(loc),'_obj',num2str(j),'_data.mat'],'firingRateSmoothing5','sumFiringRateObject5','posObjects5','amplitudeNAll','object1','sumFiringRateObject_2nd');
            saveas(gcf,[savedir1{find(unique(startlocation)==loc)},'\','neuron comparingAmplitudeNorm_binsize',num2str(binsize),'_',temp,'trial',num2str(loc),'_obj',num2str(j),'.fig'],'fig');
            saveas(gcf,[savedir1{find(unique(startlocation)==loc)},'\','neuron comparingAmplitudeNorm_binsize',num2str(binsize),'_',temp,'trial',num2str(loc),'_obj',num2str(j),'.tif'],'tif');
            saveas(gcf,[savedir1{find(unique(startlocation)==loc)},'\','neuron comparingAmplitudeNorm_binsize',num2str(binsize),'_',temp,'trial',num2str(loc),'_obj',num2str(j),'.eps'],'epsc');
        end
        close all
    end
    
%% single cell plot, before reach obj
    load([fileloc,'\','neuronIndividuals_new.mat']);
    load([fileloc,'\','further_processed_neuron_extraction_final_result.mat']);
    
    temp='S';
    
    maxS = max(neuron.C,[],2);
    thresh=3*std(neuron.S,[],2);
    binsize=5;

    for loc=unique(startlocation)
        firingrateAll={};
        countAll={};
        countTimeAll={};
        amplitudeAll={};

        counttt=1;
        for i=1:length(startlocation)
            if startlocation(i)==loc
                object1=objects{i};
                object1(:,2)=ROIlist{i}(4)-(object1(:,2));

                for j=1:size(object1,1)
                    n=Sources2D();
                    n.C=neuronIndividuals_new{i}.C(:,1:reach_obj_point_all{i}{j}/2);
                    n.S=neuronIndividuals_new{i}.S(:,1:reach_obj_point_all{i}{j}/2);
                    n.trace=neuronIndividuals_new{i}.trace(:,1:reach_obj_point_all{i}{j}/2);
                    n.time=neuronIndividuals_new{i}.time(1:reach_obj_point_all{i}{j});

                    [firingrate,count,~,countTime] = calculatingCellSpatialForSingleData_Suoqin(n,posBeforeObj{i,j},timeBeforeObj{i,j},ROIlist{i},binsize,1:size(n.C,1),thresh,temp,[],[],[0 1000000]);
                    [~,~,~,~,amplitude,amplitude_normalized] = calculatingCellSpatialForSingleData_Suoqin(n,posBeforeObj{i,j},timeBeforeObj{i,j},ROIlist{i},binsize,1:size(n.C,1),thresh,temp,[],[],[0 1000000]);  %%%bin size suggests to be 15
    %                 [firingrate,count,~,countTime] = calculatingCellSpatialForSingleData_Suoqin(neuronIndividuals_new{i},originalPos{i},originalTime{i},ROIlist{i},binsize,1:size(neuronIndividuals_new{i}.C,1),thresh,temp,[],[],[0 1000000]);
    %                 [~,~,~,~,amplitude,amplitude_normalized] = calculatingCellSpatialForSingleData_Suoqin(neuronIndividuals_new{i},originalPos{i},originalTime{i},ROIlist{i},binsize,1:size(neuronIndividuals_new{i}.C,1),thresh,temp,[],[],[0 1000000]);  %%%bin size suggests to be 15
    %                 
                    for ipp=0:length(firingrate)-1
                        firingrateAll{countt+ipp,j}=firingrate{ipp+1};
                        countAll{countt+ipp,j}=count{ipp+1};
                        countTimeAll{countt+ipp,j}=countTime;
                        amplitudeAll{countt+ipp,j}=amplitude{ipp+1};
                        amplitudeNAll{countt+ipp,j}=amplitude_normalized{ipp+1};
                    end                                         
                end
                countt=countt+length(firingrate);
            end
        end
        
        for i=1:size(countTimeAll,1)
            size1(i)=size(countTimeAll{i,1},1);
            size2(i)=size(countTimeAll{i,1},2);
        end
        size1m=max(size1);
        size2m=max(size2);
        countTime={};
        for j=1:size(countTimeAll,2)
            countTime{j}=zeros(size1m,size2m);
        end
        for i=1:size(countTimeAll,1)
            for j=1:size(countTimeAll,2)
                if ~isempty(countTimeAll{i,j})
                    countTime{j}=countTime{j}+imresize(countTimeAll{i,j},[size1m,size2m]);
                end
                if ~isempty(firingrateAll{i,j})
                    firingrateAll{i,j}=imresize(firingrateAll{i,j},[size1m,size2m]);
                end
                if ~isempty(countAll{i,j})
                    countAll{i,j}=imresize(countAll{i,j},[size1m,size2m]);
                end
                if ~isempty(amplitudeAll{i,j})
                    amplitudeAll{i,j}=imresize(amplitudeAll{i,j},[size1m,size2m]);
                end
                if ~isempty(amplitudeNAll{i,j})
                    amplitudeNAll{i,j}=imresize(amplitudeNAll{i,j},[size1m,size2m]);
                end                
            end
        end          
        
        for j=1:size(countTimeAll,2)
            plottingFiringBehaviorSpatialForSingleData_adapted_dryland({n},posBeforeObj(find(startlocation==loc),j),timeBeforeObj(find(startlocation==loc),j),ROIlist{1},countAll,[1:size(n.C,1)],thresh,[],[],savedir1(1:length(unique(startlocation))),find(unique(startlocation)==loc),5,temp,object1,firingrateAll,'firingrate',countTime,loc,j)
        end
    end

%% combine behav track plot, after reach obj
    savedirmodifier=length(unique(startlocation));
%     colorTraj = distinguishable_colors(length(behavname));
    for loc=unique(startlocation)
        figure;
        imshow(frame_bg);
        hold on;
        for i=1:length(startlocation)
            if startlocation(i)==loc
                colorTraj = distinguishable_colors(size(posBeforeObj,2));
                for j=1:size(objects{i},1)
                    plot((posAfterObj{i,j}(:,1)+behav.ROI(1))*behav.ROI3/behav.trackLength,(posAfterObj{i,j}(:,2)+behav.ROI(2))*behav.ROI3/behav.trackLength,'-','color',colorTraj(i,:));
    %                 drawnow;
    %                 pause(0.5);
                    firstLoc=min(find(~isnan(posAfterObj{i,j}(:,1))));
                    plot((posAfterObj{i,j}(firstLoc,1)+behav.ROI(1))*behav.ROI3/behav.trackLength,(posAfterObj{i,j}(firstLoc,2)+behav.ROI(2))*behav.ROI3/behav.trackLength,'c.','markerSize',20);
%                     text((posBeforeObj{i,j}(firstLoc,1)+behav.ROI(1))*behav.ROI3/behav.trackLength+2,(posBeforeObj{i,j}(firstLoc,2)+behav.ROI(2))*behav.ROI3/behav.trackLength+2,num2str(i),'color','w');
     %                 plot((originalPos{i}(:,1)+behav.ROI(1))*behav.ROI3/behav.trackLength,(originalPos{i}(:,2)+behav.ROI(2))*behav.ROI3/behav.trackLength,'-','color',colorTraj(i,:));            
                    plot((objects{i}(j,1)+behav.ROI(1))*behav.ROI3/behav.trackLength,(objects{i}(j,2)+behav.ROI(2))*behav.ROI3/behav.trackLength,'c.','markerSize',40)
                end
            end
        end
        saveas(gcf,[savedir1{find(unique(startlocation)==loc)+savedirmodifier},'\','traj_before_reach_obj_start_',num2str(loc),'.fig']);
        saveas(gcf,[savedir1{find(unique(startlocation)==loc)+savedirmodifier},'\','traj_before_reach_obj_start_',num2str(loc),'.eps'],'epsc');
        saveas(gcf,[savedir1{find(unique(startlocation)==loc)+savedirmodifier},'\','traj_before_reach_obj_start_',num2str(loc),'.tif']);
    end
    
%% combine behav track plot and neuron, after reach obj
    load([fileloc,'\','neuronIndividuals_new.mat']);
    load([fileloc,'\','further_processed_neuron_extraction_final_result.mat']);
    
    temp='S';
    
%     maxS = max(neuron.C,[],2);
    thresh =3*std(neuron.S,[],2);
    binsize=5;

    for loc=unique(startlocation)
        firingrateAll={};
        countAll={};
        countTimeAll={};
        amplitudeAll={};
        countt=1;
        for i=1:length(startlocation)
            if startlocation(i)==loc
                object1=objects{i};
                object1(:,2)=ROIlist{i}(4)-(object1(:,2));

                for j=1:size(object1,1)
                    n=Sources2D();
                    n.C=neuronIndividuals_new{i}.C(:,reach_obj_point_all{i}{j}/2:min(floor(away_from_obj_point_all{i}{j})/2-1,length(neuronIndividuals_new{i}.C)));
                    n.S=neuronIndividuals_new{i}.S(:,reach_obj_point_all{i}{j}/2:min(floor(away_from_obj_point_all{i}{j})/2-1,length(neuronIndividuals_new{i}.S)));
                    n.trace=neuronIndividuals_new{i}.trace(:,reach_obj_point_all{i}{j}/2:min(floor(away_from_obj_point_all{i}{j})/2-1,length(neuronIndividuals_new{i}.trace)));
                    n.time=neuronIndividuals_new{i}.time(reach_obj_point_all{i}{j}:min(floor(away_from_obj_point_all{i}{j})-1,length(neuronIndividuals_new{i}.time)));

                    [firingrate,count,~,countTime] = calculatingCellSpatialForSingleData_Suoqin(n,posAfterObj{i,j},timeAfterObj{i,j},ROIlist{i},binsize,1:size(n.C,1),thresh,temp,[],[],[0 1000000]);
                    [~,~,~,~,amplitude,amplitude_normalized] = calculatingCellSpatialForSingleData_Suoqin(n,posAfterObj{i,j},timeAfterObj{i,j},ROIlist{i},binsize,1:size(n.C,1),thresh,temp,[],[],[0 1000000]);  %%%bin size suggests to be 15
    %                 [firingrate,count,~,countTime] = calculatingCellSpatialForSingleData_Suoqin(neuronIndividuals_new{i},originalPos{i},originalTime{i},ROIlist{i},binsize,1:size(neuronIndividuals_new{i}.C,1),thresh,temp,[],[],[0 1000000]);
    %                 [~,~,~,~,amplitude,amplitude_normalized] = calculatingCellSpatialForSingleData_Suoqin(neuronIndividuals_new{i},originalPos{i},originalTime{i},ROIlist{i},binsize,1:size(neuronIndividuals_new{i}.C,1),thresh,temp,[],[],[0 1000000]);  %%%bin size suggests to be 15
    %                 
                    for ipp=0:length(firingrate)-1
                        firingrateAll{countt+ipp,j}=firingrate{ipp+1};
                        countAll{countt+ipp,j}=count{ipp+1};
                        countTimeAll{countt+ipp,j}=countTime;
                        amplitudeAll{countt+ipp,j}=amplitude{ipp+1};
                        amplitudeNAll{countt+ipp,j}=amplitude_normalized{ipp+1};
                    end                                         
                end
                countt=countt+length(firingrate);
            end
        end
        
        for i=1:size(countTimeAll,1)
            size1(i)=size(countTimeAll{i,1},1);
            size2(i)=size(countTimeAll{i,1},2);
        end
        size1m=max(size1);
        size2m=max(size2);
        countTime={};
        for j=1:size(countTimeAll,2)
            countTime{j}=zeros(size1m,size2m);
        end
        for i=1:size(countTimeAll,1)
            for j=1:size(countTimeAll,2)
                if ~isempty(countTimeAll{i,j})
                    countTime{j}=countTime{j}+imresize(countTimeAll{i,j},[size1m,size2m]);
                end
                if ~isempty(firingrateAll{i,j})
                    firingrateAll{i,j}=imresize(firingrateAll{i,j},[size1m,size2m]);
                end
                if ~isempty(countAll{i,j})
                    countAll{i,j}=imresize(countAll{i,j},[size1m,size2m]);
                end
                if ~isempty(amplitudeAll{i,j})
                    amplitudeAll{i,j}=imresize(amplitudeAll{i,j},[size1m,size2m]);
                end
                if ~isempty(amplitudeNAll{i,j})
                    amplitudeNAll{i,j}=imresize(amplitudeNAll{i,j},[size1m,size2m]);
                end                
            end
        end        
        
        for j=1:size(countTimeAll,2)
            [firingRateSmoothing2,sumFiringRateObject2,firingRateSmoothingt,~,radius,posObjects2,sumFiringRateObject_individual,~,~,sumFiringRateObject_2nd]=comparingFiringRateSingleConditionMultiObjects(countTimeAll(:,j)',binsize,countTime{j}, object1,'count Time',[],1,1,{'1','2','3','4','5','6','7','8','9','10'});
            save([savedir1{find(unique(startlocation)==loc)+savedirmodifier},'\','neuron_comparingCountTime_binsize',num2str(binsize),'_',temp,'trial',num2str(loc),'_obj',num2str(j),'_data.mat'],'firingRateSmoothing2','sumFiringRateObject2','posObjects2','countTimeAll','object1','sumFiringRateObject_2nd');
            saveas(gcf,[savedir1{find(unique(startlocation)==loc)+savedirmodifier},'\','neuron comparingCountTime_binsize',num2str(binsize),'_',temp,'trial',num2str(loc),'_obj',num2str(j),'.fig'],'fig');
            saveas(gcf,[savedir1{find(unique(startlocation)==loc)+savedirmodifier},'\','neuron comparingCountTime_binsize',num2str(binsize),'_',temp,'trial',num2str(loc),'_obj',num2str(j),'.tif'],'tif');
            saveas(gcf,[savedir1{find(unique(startlocation)==loc)+savedirmodifier},'\','neuron comparingCountTime_binsize',num2str(binsize),'_',temp,'trial',num2str(loc),'_obj',num2str(j),'.eps'],'epsc');
            [firingRateSmoothing3,sumFiringRateObject3,firingRateSmoothingt,~,radius,posObjects3,sumFiringRateObject_individual,~,~,sumFiringRateObject_2nd]=comparingFiringRateSingleConditionMultiObjects(firingrateAll(:,j)',binsize,countTime{j}, object1,'firing rate',[],1,1,{'1','2','3','4','5','6','7','8','9','10'});
            save([savedir1{find(unique(startlocation)==loc)+savedirmodifier},'\','neuron_comparingFiringRate_binsize',num2str(binsize),'_',temp,'trial',num2str(loc),'_obj',num2str(j),'_data.mat'],'firingRateSmoothing3','sumFiringRateObject3','posObjects3','firingrateAll','object1','sumFiringRateObject_2nd');
            saveas(gcf,[savedir1{find(unique(startlocation)==loc)+savedirmodifier},'\','neuron comparingFiringRate_binsize',num2str(binsize),'_',temp,'trial',num2str(loc),'_obj',num2str(j),'.fig'],'fig');
            saveas(gcf,[savedir1{find(unique(startlocation)==loc)+savedirmodifier},'\','neuron comparingFiringRate_binsize',num2str(binsize),'_',temp,'trial',num2str(loc),'_obj',num2str(j),'.tif'],'tif');
            saveas(gcf,[savedir1{find(unique(startlocation)==loc)+savedirmodifier},'\','neuron comparingFiringRate_binsize',num2str(binsize),'_',temp,'trial',num2str(loc),'_obj',num2str(j),'.eps'],'epsc');      
            [firingRateSmoothing5,sumFiringRateObject5,firingRateSmoothingt,~,radius,posObjects5,sumFiringRateObject_individual,~,~,sumFiringRateObject_2nd]=comparingFiringRateSingleConditionMultiObjects(amplitudeNAll(:,j)',binsize,countTime{j}, object1,'Amplitude',[],1,1,{'1','2','3','4','5','6','7','8','9','10'});
            save([savedir1{find(unique(startlocation)==loc)+savedirmodifier},'\','neuron_comparingAmplitudeNorm_binsize',num2str(binsize),'_',temp,'trial',num2str(loc),'_obj',num2str(j),'_data.mat'],'firingRateSmoothing5','sumFiringRateObject5','posObjects5','amplitudeNAll','object1','sumFiringRateObject_2nd');
            saveas(gcf,[savedir1{find(unique(startlocation)==loc)+savedirmodifier},'\','neuron comparingAmplitudeNorm_binsize',num2str(binsize),'_',temp,'trial',num2str(loc),'_obj',num2str(j),'.fig'],'fig');
            saveas(gcf,[savedir1{find(unique(startlocation)==loc)+savedirmodifier},'\','neuron comparingAmplitudeNorm_binsize',num2str(binsize),'_',temp,'trial',num2str(loc),'_obj',num2str(j),'.tif'],'tif');
            saveas(gcf,[savedir1{find(unique(startlocation)==loc)+savedirmodifier},'\','neuron comparingAmplitudeNorm_binsize',num2str(binsize),'_',temp,'trial',num2str(loc),'_obj',num2str(j),'.eps'],'epsc');
        end
        close all
    end
    
%% single cell plot, after reach obj
    load([fileloc,'\','neuronIndividuals_new.mat']);
    load([fileloc,'\','further_processed_neuron_extraction_final_result.mat']);
    
    temp='S';
    
    maxS = max(neuron.C,[],2);
    thresh=3*std(neuron.S,[],2);
    binsize=5;

    for loc=unique(startlocation)
        firingrateAll={};
        countAll={};
        countTimeAll={};
        amplitudeAll={};
        countt=1;
        for i=1:length(startlocation)
            if startlocation(i)==loc
                object1=objects{i};
                object1(:,2)=ROIlist{i}(4)-(object1(:,2));

                for j=1:size(object1,1)
                    n=Sources2D();
                    n.C=neuronIndividuals_new{i}.C(:,reach_obj_point_all{i}{j}/2:min(floor(away_from_obj_point_all{i}{j})/2-1,length(neuronIndividuals_new{i}.C)));
                    n.S=neuronIndividuals_new{i}.S(:,reach_obj_point_all{i}{j}/2:min(floor(away_from_obj_point_all{i}{j})/2-1,length(neuronIndividuals_new{i}.S)));
                    n.trace=neuronIndividuals_new{i}.trace(:,reach_obj_point_all{i}{j}/2:min(floor(away_from_obj_point_all{i}{j})/2-1,length(neuronIndividuals_new{i}.trace)));
                    n.time=neuronIndividuals_new{i}.time(reach_obj_point_all{i}{j}:min(floor(away_from_obj_point_all{i}{j})-1,length(neuronIndividuals_new{i}.time)));

                    [firingrate,count,~,countTime] = calculatingCellSpatialForSingleData_Suoqin(n,posAfterObj{i,j},timeAfterObj{i,j},ROIlist{i},binsize,1:size(n.C,1),thresh,temp,[],[],[0 1000000]);
                    [~,~,~,~,amplitude,amplitude_normalized] = calculatingCellSpatialForSingleData_Suoqin(n,posAfterObj{i,j},timeAfterObj{i,j},ROIlist{i},binsize,1:size(n.C,1),thresh,temp,[],[],[0 1000000]);  %%%bin size suggests to be 15
    %                 [firingrate,count,~,countTime] = calculatingCellSpatialForSingleData_Suoqin(neuronIndividuals_new{i},originalPos{i},originalTime{i},ROIlist{i},binsize,1:size(neuronIndividuals_new{i}.C,1),thresh,temp,[],[],[0 1000000]);
    %                 [~,~,~,~,amplitude,amplitude_normalized] = calculatingCellSpatialForSingleData_Suoqin(neuronIndividuals_new{i},originalPos{i},originalTime{i},ROIlist{i},binsize,1:size(neuronIndividuals_new{i}.C,1),thresh,temp,[],[],[0 1000000]);  %%%bin size suggests to be 15
    %                 
                    for ipp=0:length(firingrate)-1
                        firingrateAll{countt+ipp,j}=firingrate{ipp+1};
                        countAll{countt+ipp,j}=count{ipp+1};
                        countTimeAll{countt+ipp,j}=countTime;
                        amplitudeAll{countt+ipp,j}=amplitude{ipp+1};
                        amplitudeNAll{countt+ipp,j}=amplitude_normalized{ipp+1};
                    end                                         
                end
                countt=countt+length(firingrate);
            end
        end
        
        for i=1:size(countTimeAll,1)
            size1(i)=size(countTimeAll{i,1},1);
            size2(i)=size(countTimeAll{i,1},2);
        end
        size1m=max(size1);
        size2m=max(size2);
        countTime={};
        for j=1:size(countTimeAll,2)
            countTime{j}=zeros(size1m,size2m);
        end
        for i=1:size(countTimeAll,1)
            for j=1:size(countTimeAll,2)
                if ~isempty(countTimeAll{i,j})
                    countTime{j}=countTime{j}+imresize(countTimeAll{i,j},[size1m,size2m]);
                end
                if ~isempty(firingrateAll{i,j})
                    firingrateAll{i,j}=imresize(firingrateAll{i,j},[size1m,size2m]);
                end
                if ~isempty(countAll{i,j})
                    countAll{i,j}=imresize(countAll{i,j},[size1m,size2m]);
                end
                if ~isempty(amplitudeAll{i,j})
                    amplitudeAll{i,j}=imresize(amplitudeAll{i,j},[size1m,size2m]);
                end
                if ~isempty(amplitudeNAll{i,j})
                    amplitudeNAll{i,j}=imresize(amplitudeNAll{i,j},[size1m,size2m]);
                end                
            end
        end
        
        for j=1:size(countTimeAll,2)
            plottingFiringBehaviorSpatialForSingleData_adapted_dryland({n},posAfterObj(find(startlocation==loc)),timeAfterObj(find(startlocation==loc)),ROIlist{1},countAll,[1:size(n.C,1)],thresh,[],[],savedir1(1+savedirmodifier:length(unique(startlocation))+savedirmodifier),find(unique(startlocation)==loc),5,temp,object1,firingrateAll,'firingrate',countTime,loc,j)
        end
    end
    
%% combine behav track plot, all obj

    savedirmodifier=2*length(unique(startlocation));
%     colorTraj = distinguishable_colors(length(behavname));
    for loc=unique(startlocation)
        figure;
        imshow(frame_bg);
        hold on;
        colorTraj = distinguishable_colors(length(unique(startlocation)));
        for i=1:length(startlocation)
            if startlocation(i)==loc
                
                plot((originalPos{i}(:,1)+behav.ROI(1))*behav.ROI3/behav.trackLength,(originalPos{i}(:,2)+behav.ROI(2))*behav.ROI3/behav.trackLength,'-','color',colorTraj(i,:));
%                 drawnow;
%                 pause(0.5);
                firstLoc=min(find(~isnan(originalPos{i}(:,1))));
                plot((originalPos{i}(firstLoc,1)+behav.ROI(1))*behav.ROI3/behav.trackLength,(originalPos{i}(firstLoc,2)+behav.ROI(2))*behav.ROI3/behav.trackLength,'c.','markerSize',20);
%                     text((posBeforeObj{i,j}(firstLoc,1)+behav.ROI(1))*behav.ROI3/behav.trackLength+2,(posBeforeObj{i,j}(firstLoc,2)+behav.ROI(2))*behav.ROI3/behav.trackLength+2,num2str(i),'color','w');
%                 plot((originalPos{i}(:,1)+behav.ROI(1))*behav.ROI3/behav.trackLength,(originalPos{i}(:,2)+behav.ROI(2))*behav.ROI3/behav.trackLength,'-','color',colorTraj(i,:));            
                for j=1:size(objects{i},1)
                    plot((objects{i}(j,1)+behav.ROI(1))*behav.ROI3/behav.trackLength,(objects{i}(j,2)+behav.ROI(2))*behav.ROI3/behav.trackLength,'c.','markerSize',40)
                end
            end
        end
        saveas(gcf,[savedir1{find(unique(startlocation)==loc)+savedirmodifier},'\','traj_all_obj_start_',num2str(loc),'.fig']);
        saveas(gcf,[savedir1{find(unique(startlocation)==loc)+savedirmodifier},'\','traj_all_obj_start_',num2str(loc),'.eps'],'epsc');
        saveas(gcf,[savedir1{find(unique(startlocation)==loc)+savedirmodifier},'\','traj_all_obj_start_',num2str(loc),'.tif']);
    end
    
%% combine behav track plot and neuron, all obj
    load([fileloc,'\','neuronIndividuals_new.mat']);
    load([fileloc,'\','further_processed_neuron_extraction_final_result.mat']);
    
    temp='C';
    
%     maxS = max(neuron.C,[],2);
    thresh =3*std(neuron.S,[],2);
    binsize=5;

    for loc=unique(startlocation)
        firingrateAll={};
        countAll={};
        countTimeAll={};
        amplitudeAll={};
        for i=1:length(startlocation)
            if startlocation(i)==loc
                object1=objects{i};
                object1(:,2)=ROIlist{i}(4)-(object1(:,2));
                n=Sources2D();
                n.C=neuronIndividuals_new{i}.C;
                n.S=neuronIndividuals_new{i}.S;
                n.trace=neuronIndividuals_new{i}.trace;
                n.time=neuronIndividuals_new{i}.time;

                [firingrate,count,~,countTime] = calculatingCellSpatialForSingleData_Suoqin(n,originalPos{i},originalTime{i},ROIlist{i},binsize,1:size(n.C,1),thresh,temp,[],[],[0 1000000]);
                [~,~,~,~,amplitude,amplitude_normalized] = calculatingCellSpatialForSingleData_Suoqin(n,originalPos{i},originalTime{i},ROIlist{i},binsize,1:size(n.C,1),thresh,temp,[],[],[0 1000000]);  %%%bin size suggests to be 15
%                 
                for ipp=1:length(firingrate)
                    firingrateAll{ipp}=firingrate{ipp};
                    countAll{ipp}=count{ipp};
                    countTimeAll{ipp}=countTime;
                    amplitudeAll{ipp}=amplitude{ipp};
                    amplitudeNAll{ipp}=amplitude_normalized{ipp};
                end                                         
            end
        end
        
        for i=1:size(countTimeAll,1)
            size1(i)=size(countTimeAll{i,1},1);
            size2(i)=size(countTimeAll{i,1},2);
        end
        size1m=max(size1);
        size2m=max(size2);
        countTime=zeros(size1m,size2m);

        for i=1:size(countTimeAll,1)
            if ~isempty(countTimeAll{i})
                countTime=countTime+imresize(countTimeAll{i},[size1m,size2m]);
            end
            if ~isempty(firingrateAll{i})
                firingrateAll{i}=imresize(firingrateAll{i},[size1m,size2m]);
            end
            if ~isempty(countAll{i})
                countAll{i}=imresize(countAll{i},[size1m,size2m]);
            end
            if ~isempty(amplitudeAll{i})
                amplitudeAll{i}=imresize(amplitudeAll{i},[size1m,size2m]);
            end
            if ~isempty(amplitudeNAll{i})
                amplitudeNAll{i}=imresize(amplitudeNAll{i},[size1m,size2m]);
            end                
        end        

        [firingRateSmoothing2,sumFiringRateObject2,firingRateSmoothingt,~,radius,posObjects2,sumFiringRateObject_individual,~,~,sumFiringRateObject_2nd]=comparingFiringRateSingleConditionMultiObjects(countTimeAll,binsize,countTime, object1,'count Time',[],1,1,{'1','2','3','4','5','6','7','8','9','10'});
        save([savedir1{find(unique(startlocation)==loc)+savedirmodifier},'\','neuron_comparingCountTime_binsize',num2str(binsize),'_',temp,'trial',num2str(loc),'_obj',num2str(j),'_data.mat'],'firingRateSmoothing2','sumFiringRateObject2','posObjects2','countTimeAll','object1','sumFiringRateObject_2nd');
        saveas(gcf,[savedir1{find(unique(startlocation)==loc)+savedirmodifier},'\','neuron comparingCountTime_binsize',num2str(binsize),'_',temp,'trial',num2str(loc),'_obj',num2str(j),'.fig'],'fig');
        saveas(gcf,[savedir1{find(unique(startlocation)==loc)+savedirmodifier},'\','neuron comparingCountTime_binsize',num2str(binsize),'_',temp,'trial',num2str(loc),'_obj',num2str(j),'.tif'],'tif');
        saveas(gcf,[savedir1{find(unique(startlocation)==loc)+savedirmodifier},'\','neuron comparingCountTime_binsize',num2str(binsize),'_',temp,'trial',num2str(loc),'_obj',num2str(j),'.eps'],'epsc');
        [firingRateSmoothing3,sumFiringRateObject3,firingRateSmoothingt,~,radius,posObjects3,sumFiringRateObject_individual,~,~,sumFiringRateObject_2nd]=comparingFiringRateSingleConditionMultiObjects(firingrateAll,binsize,countTime, object1,'firing rate',[],1,1,{'1','2','3','4','5','6','7','8','9','10'});
        save([savedir1{find(unique(startlocation)==loc)+savedirmodifier},'\','neuron_comparingFiringRate_binsize',num2str(binsize),'_',temp,'trial',num2str(loc),'_obj',num2str(j),'_data.mat'],'firingRateSmoothing3','sumFiringRateObject3','posObjects3','firingrateAll','object1','sumFiringRateObject_2nd');
        saveas(gcf,[savedir1{find(unique(startlocation)==loc)+savedirmodifier},'\','neuron comparingFiringRate_binsize',num2str(binsize),'_',temp,'trial',num2str(loc),'_obj',num2str(j),'.fig'],'fig');
        saveas(gcf,[savedir1{find(unique(startlocation)==loc)+savedirmodifier},'\','neuron comparingFiringRate_binsize',num2str(binsize),'_',temp,'trial',num2str(loc),'_obj',num2str(j),'.tif'],'tif');
        saveas(gcf,[savedir1{find(unique(startlocation)==loc)+savedirmodifier},'\','neuron comparingFiringRate_binsize',num2str(binsize),'_',temp,'trial',num2str(loc),'_obj',num2str(j),'.eps'],'epsc');      
        [firingRateSmoothing5,sumFiringRateObject5,firingRateSmoothingt,~,radius,posObjects5,sumFiringRateObject_individual,~,~,sumFiringRateObject_2nd]=comparingFiringRateSingleConditionMultiObjects(amplitudeNAll,binsize,countTime, object1,'Amplitude',[],1,1,{'1','2','3','4','5','6','7','8','9','10'});
        save([savedir1{find(unique(startlocation)==loc)+savedirmodifier},'\','neuron_comparingAmplitudeNorm_binsize',num2str(binsize),'_',temp,'trial',num2str(loc),'_obj',num2str(j),'_data.mat'],'firingRateSmoothing5','sumFiringRateObject5','posObjects5','amplitudeNAll','object1','sumFiringRateObject_2nd');
        saveas(gcf,[savedir1{find(unique(startlocation)==loc)+savedirmodifier},'\','neuron comparingAmplitudeNorm_binsize',num2str(binsize),'_',temp,'trial',num2str(loc),'_obj',num2str(j),'.fig'],'fig');
        saveas(gcf,[savedir1{find(unique(startlocation)==loc)+savedirmodifier},'\','neuron comparingAmplitudeNorm_binsize',num2str(binsize),'_',temp,'trial',num2str(loc),'_obj',num2str(j),'.tif'],'tif');
        saveas(gcf,[savedir1{find(unique(startlocation)==loc)+savedirmodifier},'\','neuron comparingAmplitudeNorm_binsize',num2str(binsize),'_',temp,'trial',num2str(loc),'_obj',num2str(j),'.eps'],'epsc');
        close all
    end
    