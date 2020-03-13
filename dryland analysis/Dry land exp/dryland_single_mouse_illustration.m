function dryland_single_mouse_illustration(foldernamet,behavname,startlocation,savedir,fileloc)

cd(foldernamet);

savedirall{1}=[savedir,'\','start_right_top_before_obj'];
savedirall{2}=[savedir,'\','start_right_bottom_before_obj'];
savedirall{3}=[savedir,'\','start_left_bottom_before_obj'];
savedirall{4}=[savedir,'\','start_left_top_before_obj'];
savedirall{5}=[savedir,'\','start_right_top_after_obj'];
savedirall{6}=[savedir,'\','start_right_bottom_after_obj'];
savedirall{7}=[savedir,'\','start_left_bottom_after_obj'];
savedirall{8}=[savedir,'\','start_left_top_after_obj'];
savedirall{9}=[savedir,'\','reward_before_obj'];
savedirall{10}=[savedir,'\','unexpected_no_reward_before_obj'];
savedirall{11}=[savedir,'\','reward_after_obj'];
savedirall{12}=[savedir,'\','unexpected_no_reward_after_obj'];
savedirall{13}=[savedir,'\','reward_before_obj_xx'];
savedirall{14}=[savedir,'\','unexpected_no_reward_before_obj_xx'];
savedirall{15}=[savedir,'\','other_reward_before_obj_xx'];
savedirall{16}=[savedir,'\','reward_after_obj_xx'];
savedirall{17}=[savedir,'\','unexpected_no_reward_after_obj_xx'];
savedirall{18}=[savedir,'\','other_reward_after_obj_xx'];

% savedir1=savedirall([13 14 15 16 17 18]);
savedir1=savedirall([1 3 4 5 7 8]);

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

objmodifier=[4 -4];%amend position

for i=1:length(behavname)
    load(behavname{i});
    object=behav.object;
    object(1)=object(1)+objmodifier(1);
    object(2)=behav.ROI(4)-(object(2)+objmodifier(2));
    behavpos=behav.position;
    behavtime=behav.time;
    coor_diff=(behavpos-repmat(object,size(behavpos,1),1));
    dis_to_obj=(coor_diff(:,1).^2+coor_diff(:,2).^2).^0.5;
    dis_to_obj(isnan(dis_to_obj))=inf;
%     reach_obj_point=find(dis_to_obj==min(dis_to_obj));
    reach_obj_point=min(find(dis_to_obj<10));
    if ~isempty(min(find(dis_to_obj(reach_obj_point:end)>40)))
        away_from_obj_point=min(find(dis_to_obj(reach_obj_point:end)>40))+reach_obj_point;
    else
        away_from_obj_point=length(behavtime);
    end

    originalPos{i}=behavpos;
    originalTime{i}=behavtime;    
    posBeforeObj{i}=behavpos(1:reach_obj_point,:);
    timeBeforeObj{i}=behavtime(1:reach_obj_point);
    posAfterObj{i}=behavpos(reach_obj_point:away_from_obj_point,:);
    timeAfterObj{i}=behavtime(reach_obj_point:away_from_obj_point);
    dis_to_obj_all{i}=dis_to_obj;
    reach_obj_point_all{i}=reach_obj_point;
    away_from_obj_point_all{i}=away_from_obj_point;
    ROIlist{i}=behav.ROI;
end

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
    colorTraj = distinguishable_colors(length(behavname));

    for loc=unique(startlocation)
        figure;
        imshow(frame_bg);
        hold on;
        for i=1:length(startlocation)
            if startlocation(i)==loc
                plot((posBeforeObj{i}(:,1)+behav.ROI(1))*behav.ROI3/behav.trackLength,(posBeforeObj{i}(:,2)+behav.ROI(2))*behav.ROI3/behav.trackLength,'-','color',colorTraj(i,:));
%                 drawnow;
%                 pause(0.5);
                firstLoc=min(find(~isnan(posBeforeObj{i}(:,1))));
                plot((posBeforeObj{i}(firstLoc,1)+behav.ROI(1))*behav.ROI3/behav.trackLength,(posBeforeObj{i}(firstLoc,2)+behav.ROI(2))*behav.ROI3/behav.trackLength,'c.','markerSize',20);
                text((posBeforeObj{i}(firstLoc,1)+behav.ROI(1))*behav.ROI3/behav.trackLength+2,(posBeforeObj{i}(firstLoc,2)+behav.ROI(2))*behav.ROI3/behav.trackLength+2,num2str(i),'color','w');
 %                 plot((originalPos{i}(:,1)+behav.ROI(1))*behav.ROI3/behav.trackLength,(originalPos{i}(:,2)+behav.ROI(2))*behav.ROI3/behav.trackLength,'-','color',colorTraj(i,:));            
            end
        end
        plot((object(1)+behav.ROI(1))*behav.ROI3/behav.trackLength,(object(2)+behav.ROI(2))*behav.ROI3/behav.trackLength,'c.','markerSize',40)
        saveas(gcf,[savedir1{find(unique(startlocation)==loc)},'\','traj_before_reach_obj_start_',num2str(loc),'.fig']);
        saveas(gcf,[savedir1{find(unique(startlocation)==loc)},'\','traj_before_reach_obj_start_',num2str(loc),'.eps'],'epsc');
        saveas(gcf,[savedir1{find(unique(startlocation)==loc)},'\','traj_before_reach_obj_start_',num2str(loc),'.tif']);
    end
    
%% combine behav track plot and neuron, before reach obj
    load([fileloc,'\','neuronIndividuals_new.mat']);
    load([fileloc,'\','further_processed_neuron_extraction_final_result.mat']);
    
    temp='C';
    for i=1:length(startlocation)
        neuronIndividuals_new{i}.C=neuronIndividuals_new{i}.C(:,1:reach_obj_point_all{i}/2);
        neuronIndividuals_new{i}.S=neuronIndividuals_new{i}.S(:,1:reach_obj_point_all{i}/2);
        neuronIndividuals_new{i}.trace=neuronIndividuals_new{i}.trace(:,1:reach_obj_point_all{i}/2);
        neuronIndividuals_new{i}.time=neuronIndividuals_new{i}.time(1:reach_obj_point_all{i});
    end
    
%     maxS = max(neuron.C,[],2);
    thresh =3*std(neuron.S,[],2);
    binsize=5;

    for loc=unique(startlocation)
        firingrateAll={};
        countAll={};
        countTimeAll={};
        amplitudeAll={};
        countt=1;
        object1=object;
        object1(2)=ROIlist{1}(4)-(object1(2));
        for i=1:length(startlocation)
            if startlocation(i)==loc
                [firingrate,count,~,countTime] = calculatingCellSpatialForSingleData_Suoqin(neuronIndividuals_new{i},posBeforeObj{i},timeBeforeObj{i},ROIlist{i},binsize,1:size(neuronIndividuals_new{i}.C,1),thresh,temp,[],[],[0 1000000]);
                [~,~,~,~,amplitude,amplitude_normalized] = calculatingCellSpatialForSingleData_Suoqin(neuronIndividuals_new{i},posBeforeObj{i},timeBeforeObj{i},ROIlist{i},binsize,1:size(neuronIndividuals_new{i}.C,1),thresh,temp,[],[],[0 1000000]);  %%%bin size suggests to be 15
%                 [firingrate,count,~,countTime] = calculatingCellSpatialForSingleData_Suoqin(neuronIndividuals_new{i},originalPos{i},originalTime{i},ROIlist{i},binsize,1:size(neuronIndividuals_new{i}.C,1),thresh,temp,[],[],[0 1000000]);
%                 [~,~,~,~,amplitude,amplitude_normalized] = calculatingCellSpatialForSingleData_Suoqin(neuronIndividuals_new{i},originalPos{i},originalTime{i},ROIlist{i},binsize,1:size(neuronIndividuals_new{i}.C,1),thresh,temp,[],[],[0 1000000]);  %%%bin size suggests to be 15
%                 
                firingrateAll(countt:countt+length(firingrate)-1)=firingrate;
                countAll(countt:countt+length(firingrate)-1)=count;
                for ipp=0:length(firingrate)-1
                    countTimeAll{countt+ipp}=countTime;
                end
                amplitudeAll(countt:countt+length(firingrate)-1)=amplitude;
                countt=countt+length(firingrate);
            end
        end
        
        for i=1:length(countTimeAll)
            size1=size(countTimeAll{i},1);
            size2=size(countTimeAll{i},2);
        end
        size1m=max(size1);
        size2m=max(size2);
        countTime=zeros(size1m,size2m);
        for i=1:length(countTimeAll)
            countTime=countTime+imresize(countTimeAll{i},[size1m,size2m]);
            if ~isempty(firingrateAll{i})
                firingrateAll{i}=imresize(firingrateAll{i},[size1m,size2m]);
            end
            if ~isempty(countAll{i})
                countAll{i}=imresize(countAll{i},[size1m,size2m]);
            end
            if ~isempty(amplitudeAll{i})
                amplitudeAll{i}=imresize(amplitudeAll{i},[size1m,size2m]);
            end
        end        
        
        [firingRateSmoothing1,sumFiringRateObject1,firingRateSmoothingt,~,radius,posObjects1,sumFiringRateObject_individual,~,~,sumFiringRateObject_2nd]=comparingFiringRateSingleConditionMultiObjects(countAll,binsize,countTime, object1,'events',[],1,1,{'obj'});
        save([savedir1{find(unique(startlocation)==loc)},'\','neuron_comparingCountevents_binsize',num2str(binsize),'_',temp,'start',num2str(loc),'.mat'],'firingRateSmoothing1','sumFiringRateObject1','posObjects1','countAll','object1','sumFiringRateObject_2nd');
        saveas(gcf,[savedir1{find(unique(startlocation)==loc)},'\','neuron comparingCountevents_binsize',num2str(binsize),'_',temp,'start',num2str(loc),'.fig'],'fig');
        saveas(gcf,[savedir1{find(unique(startlocation)==loc)},'\','neuron comparingCountevents_binsize',num2str(binsize),'_',temp,'start',num2str(loc),'.tif'],'tif');
        saveas(gcf,[savedir1{find(unique(startlocation)==loc)},'\','neuron comparingCountevents_binsize',num2str(binsize),'_',temp,'start',num2str(loc),'.eps'],'epsc');
        [firingRateSmoothing2,sumFiringRateObject2,firingRateSmoothingt,~,radius,posObjects2,sumFiringRateObject_individual,~,~,sumFiringRateObject_2nd]=comparingFiringRateSingleConditionMultiObjects(countTimeAll,binsize,countTime, object1,'count Time',[],1,1,{'obj'});
        save([savedir1{find(unique(startlocation)==loc)},'\','neuron_comparingCountTime_binsize',num2str(binsize),'_',temp,'start',num2str(loc),'_data.mat'],'firingRateSmoothing2','sumFiringRateObject2','posObjects2','countTimeAll','object1','sumFiringRateObject_2nd');
        saveas(gcf,[savedir1{find(unique(startlocation)==loc)},'\','neuron comparingCountTime_binsize',num2str(binsize),'_',temp,'start',num2str(loc),'.fig'],'fig');
        saveas(gcf,[savedir1{find(unique(startlocation)==loc)},'\','neuron comparingCountTime_binsize',num2str(binsize),'_',temp,'start',num2str(loc),'.tif'],'tif');
        saveas(gcf,[savedir1{find(unique(startlocation)==loc)},'\','neuron comparingCountTime_binsize',num2str(binsize),'_',temp,'start',num2str(loc),'.eps'],'epsc');
        [firingRateSmoothing3,sumFiringRateObject3,firingRateSmoothingt,~,radius,posObjects3,sumFiringRateObject_individual,~,~,sumFiringRateObject_2nd]=comparingFiringRateSingleConditionMultiObjects(firingrateAll,binsize,countTime, object1,'firing rate',[],1,1,{'obj'});
        save([savedir1{find(unique(startlocation)==loc)},'\','neuron_comparingFiringRate_binsize',num2str(binsize),'_',temp,'start',num2str(loc),'_data.mat'],'firingRateSmoothing3','sumFiringRateObject3','posObjects3','firingrateAll','object1','sumFiringRateObject_2nd');
        saveas(gcf,[savedir1{find(unique(startlocation)==loc)},'\','neuron comparingFiringRate_binsize',num2str(binsize),'_',temp,'start',num2str(loc),'.fig'],'fig');
        saveas(gcf,[savedir1{find(unique(startlocation)==loc)},'\','neuron comparingFiringRate_binsize',num2str(binsize),'_',temp,'start',num2str(loc),'.tif'],'tif');
        saveas(gcf,[savedir1{find(unique(startlocation)==loc)},'\','neuron comparingFiringRate_binsize',num2str(binsize),'_',temp,'start',num2str(loc),'.eps'],'epsc');      
        [firingRateSmoothing4,sumFiringRateObject4,firingRateSmoothingt,~,radius,posObjects4,sumFiringRateObject_individual,~,~,sumFiringRateObject_2nd]=comparingFiringRateSingleConditionMultiObjects(amplitudeAll,binsize,countTime, object1,'Amplitude',[],1,1,{'obj'});
        save([savedir1{find(unique(startlocation)==loc)},'\','neuron_comparingAmplitude_binsize',num2str(binsize),'_',temp,'start',num2str(loc),'_data.mat'],'firingRateSmoothing4','sumFiringRateObject4','posObjects4','amplitudeAll','object1','sumFiringRateObject_2nd');
        saveas(gcf,[savedir1{find(unique(startlocation)==loc)},'\','neuron comparingAmplitude_binsize',num2str(binsize),'_',temp,'start',num2str(loc),'.fig'],'fig');
        saveas(gcf,[savedir1{find(unique(startlocation)==loc)},'\','neuron comparingAmplitude_binsize',num2str(binsize),'_',temp,'start',num2str(loc),'.tif'],'tif');
        saveas(gcf,[savedir1{find(unique(startlocation)==loc)},'\','neuron comparingAmplitude_binsize',num2str(binsize),'_',temp,'start',num2str(loc),'.eps'],'epsc');
    end
    
%% single cell plot, before reach obj
    load([fileloc,'\','neuronIndividuals_new.mat']);
    load([fileloc,'\','further_processed_neuron_extraction_final_result.mat']);
    
    temp='C';
    for i=1:length(startlocation)
        neuronIndividuals_new{i}.C=neuronIndividuals_new{i}.C(:,1:reach_obj_point_all{i}/2);
        neuronIndividuals_new{i}.S=neuronIndividuals_new{i}.S(:,1:reach_obj_point_all{i}/2);
        neuronIndividuals_new{i}.trace=neuronIndividuals_new{i}.trace(:,1:reach_obj_point_all{i}/2);
        neuronIndividuals_new{i}.time=neuronIndividuals_new{i}.time(1:reach_obj_point_all{i});
    end
    
    maxS = max(neuron.C,[],2);
    thresh=3*std(neuron.S,[],2);
    binsize=5;

    for loc=unique(startlocation)
        firingrateAll={};
        countAll={};
        countTimeAll={};
        amplitudeAll={};
        object1=object;
        object1(2)=ROIlist{1}(4)-(object1(2));

        counttt=1;
        for i=1:length(startlocation)
            if startlocation(i)==loc
                [firingrate,count,~,countTime] = calculatingCellSpatialForSingleData_Suoqin(neuronIndividuals_new{i},posBeforeObj{i},timeBeforeObj{i},ROIlist{i},binsize,1:size(neuronIndividuals_new{i}.C,1),thresh,temp,[],[],[0 1000000]);
                [~,~,~,~,amplitude,amplitude_normalized] = calculatingCellSpatialForSingleData_Suoqin(neuronIndividuals_new{i},posBeforeObj{i},timeBeforeObj{i},ROIlist{i},binsize,1:size(neuronIndividuals_new{i}.C,1),thresh,temp,[],[],[0 1000000]);  %%%bin size suggests to be 15
%                 [firingrate,count,~,countTime] = calculatingCellSpatialForSingleData_Suoqin(neuronIndividuals_new{i},originalPos{i},originalTime{i},ROIlist{i},binsize,1:size(neuronIndividuals_new{i}.C,1),thresh,temp,[],[],[0 1000000]);
%                 [~,~,~,~,amplitude,amplitude_normalized] = calculatingCellSpatialForSingleData_Suoqin(neuronIndividuals_new{i},originalPos{i},originalTime{i},ROIlist{i},binsize,1:size(neuronIndividuals_new{i}.C,1),thresh,temp,[],[],[0 1000000]);  %%%bin size suggests to be 15
%                 
                firingrateAll(1:size(firingrate,2),counttt)=firingrate;
                countAll(1:size(firingrate,2),counttt)=count;
                for ipp=1:length(firingrate)
                    countTimeAll{ipp,counttt}=countTime;
                end
                amplitudeAll(1:size(firingrate,2),counttt)=amplitude;
                counttt=counttt+1;
            end
        end
        
        for i=1:length(countTimeAll)
            size1=size(countTimeAll{i},1);
            size2=size(countTimeAll{i},2);
        end
        size1m=max(size1);
        size2m=max(size2);
        countTime=zeros(size1m,size2m);
        for i=1:length(countTimeAll)
            for j=1:size(firingrateAll,2)
                countTime=countTime+imresize(countTimeAll{i,j},[size1m,size2m]);
                if ~isempty(firingrateAll{i,j})
                    firingrateAll{i,j}=imresize(firingrateAll{i,j},[size1m,size2m]);
                else
                    firingrateAll{i,j}=zeros(size1m,size2m);
                end
                if ~isempty(countAll{i,j})
                    countAll{i,j}=imresize(countAll{i,j},[size1m,size2m]);
                else
                    countAll{i,j}=zeros(size1m,size2m);
                end
                if ~isempty(amplitudeAll{i,j})
                    amplitudeAll{i,j}=imresize(amplitudeAll{i,j},[size1m,size2m]);
                else
                    amplitudeAll{i,j}=zeros(size1m,size2m);
                end
            end
        end     
        
        firingrateAll_individual=firingrateAll(:,1);
        countAll_individual=countAll(:,1);
        amplitudeAll_individual=amplitudeAll(:,1);
        for j=2:size(firingrateAll,2)
            for i=1:length(countTimeAll)
                if ~isempty(firingrateAll{i,j})
                    firingrateAll_individual{i}=firingrateAll_individual{i}+firingrateAll{i,j};
                end
                if ~isempty(countAll{i,j})
                    countAll_individual{i}=countAll_individual{i}+countAll{i,j};
                end
                if ~isempty(amplitudeAll{i,j})
                    amplitudeAll_individual{i}=amplitudeAll_individual{i}+amplitudeAll{i,j};
                end
            end
        end
        
        plottingFiringBehaviorSpatialForSingleData_adapted_dryland(neuronIndividuals_new(find(startlocation==loc)),posBeforeObj(find(startlocation==loc)),timeBeforeObj(find(startlocation==loc)),ROIlist{1},countAll_individual,[1:size(neuron.C,1)],thresh,[],[],savedir1(1:3),find(unique(startlocation)==loc),5,'C',object1,firingrateAll_individual,'firingrate',countTime,loc)
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
                plot((posAfterObj{i}(:,1)+behav.ROI(1))*behav.ROI3/behav.trackLength,(posAfterObj{i}(:,2)+behav.ROI(2))*behav.ROI3/behav.trackLength,'-','color',colorTraj(i,:));
%                 drawnow;
%                 pause(0.5);
                firstLoc=find(~isnan(posAfterObj{i}(:,1)), 1 );
                plot((posAfterObj{i}(firstLoc,1)+behav.ROI(1))*behav.ROI3/behav.trackLength,(posAfterObj{i}(firstLoc,2)+behav.ROI(2))*behav.ROI3/behav.trackLength,'c.','markerSize',20);
%                 text((posBeforeObj{i}(firstLoc,1)+behav.ROI(1))*behav.ROI3/behav.trackLength,(posBeforeObj{i}(firstLoc,2)+behav.ROI(2))*behav.ROI3/behav.trackLength,num2str(i),'color','w');
%                 plot((originalPos{i}(:,1)+behav.ROI(1))*behav.ROI3/behav.trackLength,(originalPos{i}(:,2)+behav.ROI(2))*behav.ROI3/behav.trackLength,'-','color',colorTraj(i,:));            
            end
        end
        plot((object(1)+behav.ROI(1))*behav.ROI3/behav.trackLength,(object(2)+behav.ROI(2))*behav.ROI3/behav.trackLength,'c.','markerSize',40)
        saveas(gcf,[savedir1{find(unique(startlocation)==loc)+savedirmodifier},'\','traj_After_reach_obj_start_',num2str(loc),'.fig']);
        saveas(gcf,[savedir1{find(unique(startlocation)==loc)+savedirmodifier},'\','traj_After_reach_obj_start_',num2str(loc),'.eps'],'epsc');
        saveas(gcf,[savedir1{find(unique(startlocation)==loc)+savedirmodifier},'\','traj_After_reach_obj_start_',num2str(loc),'.tif']);
    end
    
%% combine behav track plot and neuron, after reach obj
    load([fileloc,'\','neuronIndividuals_new.mat']);
    load([fileloc,'\','further_processed_neuron_extraction_final_result.mat']);
    
    temp='C';
    for i=1:length(startlocation)
        neuronIndividuals_new{i}.C=neuronIndividuals_new{i}.C(:,reach_obj_point_all{i}/2:min(floor(away_from_obj_point_all{i})/2-1,length(neuronIndividuals_new{i}.C)));
        neuronIndividuals_new{i}.S=neuronIndividuals_new{i}.S(:,reach_obj_point_all{i}/2:min(floor(away_from_obj_point_all{i})/2-1,length(neuronIndividuals_new{i}.S)));
        neuronIndividuals_new{i}.trace=neuronIndividuals_new{i}.trace(:,reach_obj_point_all{i}/2:min(floor(away_from_obj_point_all{i})/2-1,length(neuronIndividuals_new{i}.trace)));
        neuronIndividuals_new{i}.time=neuronIndividuals_new{i}.time(reach_obj_point_all{i}:min(floor(away_from_obj_point_all{i})-1,length(neuronIndividuals_new{i}.time)));
    end
    
    maxS = max(neuron.C,[],2);
    thresh=3*std(neuron.S,[],2);
    binsize=5;

    for loc=unique(startlocation)
        firingrateAll={};
        countAll={};
        countTimeAll={};
        amplitudeAll={};
        countt=1;
        object1=object;
        object1(2)=ROIlist{1}(4)-(object1(2));
        for i=1:length(startlocation)
            if startlocation(i)==loc
                [firingrate,count,~,countTime] = calculatingCellSpatialForSingleData_Suoqin(neuronIndividuals_new{i},posAfterObj{i},timeAfterObj{i},ROIlist{i},binsize,1:size(neuronIndividuals_new{i}.C,1),thresh,temp,[],[],[0 1000000]);
                [~,~,~,~,amplitude,amplitude_normalized] = calculatingCellSpatialForSingleData_Suoqin(neuronIndividuals_new{i},posAfterObj{i},timeAfterObj{i},ROIlist{i},binsize,1:size(neuronIndividuals_new{i}.C,1),thresh,temp,[],[],[0 1000000]);  %%%bin size suggests to be 15
%                 [firingrate,count,~,countTime] = calculatingCellSpatialForSingleData_Suoqin(neuronIndividuals_new{i},originalPos{i},originalTime{i},ROIlist{i},binsize,1:size(neuronIndividuals_new{i}.C,1),thresh,temp,[],[],[0 1000000]);
%                 [~,~,~,~,amplitude,amplitude_normalized] = calculatingCellSpatialForSingleData_Suoqin(neuronIndividuals_new{i},originalPos{i},originalTime{i},ROIlist{i},binsize,1:size(neuronIndividuals_new{i}.C,1),thresh,temp,[],[],[0 1000000]);  %%%bin size suggests to be 15
%                 
                firingrateAll(countt:countt+length(firingrate)-1)=firingrate;
                countAll(countt:countt+length(firingrate)-1)=count;
                for ipp=0:length(firingrate)-1
                    countTimeAll{countt+ipp}=countTime;
                end
                amplitudeAll(countt:countt+length(firingrate)-1)=amplitude;
                countt=countt+length(firingrate);
            end
        end
        
        for i=1:length(countTimeAll)
            size1=size(countTimeAll{i},1);
            size2=size(countTimeAll{i},2);
        end
        size1m=max(size1);
        size2m=max(size2);
        countTime=zeros(size1m,size2m);
        for i=1:length(countTimeAll)
            countTime=countTime+imresize(countTimeAll{i},[size1m,size2m]);
            if ~isempty(firingrateAll{i})
                firingrateAll{i}=imresize(firingrateAll{i},[size1m,size2m]);
            end
            if ~isempty(countAll{i})
                countAll{i}=imresize(countAll{i},[size1m,size2m]);
            end
            if ~isempty(amplitudeAll{i})
                amplitudeAll{i}=imresize(amplitudeAll{i},[size1m,size2m]);
            end
        end        
        
        [firingRateSmoothing1,sumFiringRateObject1,firingRateSmoothingt,~,radius,posObjects1,sumFiringRateObject_individual,~,~,sumFiringRateObject_2nd]=comparingFiringRateSingleConditionMultiObjects(countAll,binsize,countTime, object1,'events',[],1,1,{'obj'});
        save([savedir1{find(unique(startlocation)==loc)+savedirmodifier},'\','neuron_comparingCountevents_binsize',num2str(binsize),'_',temp,'start',num2str(loc),'.mat'],'firingRateSmoothing1','sumFiringRateObject1','posObjects1','countAll','object1','sumFiringRateObject_2nd');
        saveas(gcf,[savedir1{find(unique(startlocation)==loc)+savedirmodifier},'\','neuron comparingCountevents_binsize',num2str(binsize),'_',temp,'start',num2str(loc),'.fig'],'fig');
        saveas(gcf,[savedir1{find(unique(startlocation)==loc)+savedirmodifier},'\','neuron comparingCountevents_binsize',num2str(binsize),'_',temp,'start',num2str(loc),'.tif'],'tif');
        saveas(gcf,[savedir1{find(unique(startlocation)==loc)+savedirmodifier},'\','neuron comparingCountevents_binsize',num2str(binsize),'_',temp,'start',num2str(loc),'.eps'],'epsc');
        [firingRateSmoothing2,sumFiringRateObject2,firingRateSmoothingt,~,radius,posObjects2,sumFiringRateObject_individual,~,~,sumFiringRateObject_2nd]=comparingFiringRateSingleConditionMultiObjects(countTimeAll,binsize,countTime, object1,'count Time',[],1,1,{'obj'});
        save([savedir1{find(unique(startlocation)==loc)+savedirmodifier},'\','neuron_comparingCountTime_binsize',num2str(binsize),'_',temp,'start',num2str(loc),'_data.mat'],'firingRateSmoothing2','sumFiringRateObject2','posObjects2','countTimeAll','object1','sumFiringRateObject_2nd');
        saveas(gcf,[savedir1{find(unique(startlocation)==loc)+savedirmodifier},'\','neuron comparingCountTime_binsize',num2str(binsize),'_',temp,'start',num2str(loc),'.fig'],'fig');
        saveas(gcf,[savedir1{find(unique(startlocation)==loc)+savedirmodifier},'\','neuron comparingCountTime_binsize',num2str(binsize),'_',temp,'start',num2str(loc),'.tif'],'tif');
        saveas(gcf,[savedir1{find(unique(startlocation)==loc)+savedirmodifier},'\','neuron comparingCountTime_binsize',num2str(binsize),'_',temp,'start',num2str(loc),'.eps'],'epsc');
        [firingRateSmoothing3,sumFiringRateObject3,firingRateSmoothingt,~,radius,posObjects3,sumFiringRateObject_individual,~,~,sumFiringRateObject_2nd]=comparingFiringRateSingleConditionMultiObjects(firingrateAll,binsize,countTime, object1,'firing rate',[],1,1,{'obj'});
        save([savedir1{find(unique(startlocation)==loc)},'\','neuron_comparingFiringRate_binsize',num2str(binsize),'_',temp,'start',num2str(loc),'_data.mat'],'firingRateSmoothing3','sumFiringRateObject3','posObjects3','firingrateAll','object1','sumFiringRateObject_2nd');
        saveas(gcf,[savedir1{find(unique(startlocation)==loc)+savedirmodifier},'\','neuron comparingFiringRate_binsize',num2str(binsize),'_',temp,'start',num2str(loc),'.fig'],'fig');
        saveas(gcf,[savedir1{find(unique(startlocation)==loc)+savedirmodifier},'\','neuron comparingFiringRate_binsize',num2str(binsize),'_',temp,'start',num2str(loc),'.tif'],'tif');
        saveas(gcf,[savedir1{find(unique(startlocation)==loc)+savedirmodifier},'\','neuron comparingFiringRate_binsize',num2str(binsize),'_',temp,'start',num2str(loc),'.eps'],'epsc');      
        [firingRateSmoothing4,sumFiringRateObject4,firingRateSmoothingt,~,radius,posObjects4,sumFiringRateObject_individual,~,~,sumFiringRateObject_2nd]=comparingFiringRateSingleConditionMultiObjects(amplitudeAll,binsize,countTime, object1,'Amplitude',[],1,1,{'obj'});
        save([savedir1{find(unique(startlocation)==loc)+savedirmodifier},'\','neuron_comparingAmplitude_binsize',num2str(binsize),'_',temp,'start',num2str(loc),'_data.mat'],'firingRateSmoothing4','sumFiringRateObject4','posObjects4','amplitudeAll','object1','sumFiringRateObject_2nd');
        saveas(gcf,[savedir1{find(unique(startlocation)==loc)+savedirmodifier},'\','neuron comparingAmplitude_binsize',num2str(binsize),'_',temp,'start',num2str(loc),'.fig'],'fig');
        saveas(gcf,[savedir1{find(unique(startlocation)==loc)+savedirmodifier},'\','neuron comparingAmplitude_binsize',num2str(binsize),'_',temp,'start',num2str(loc),'.tif'],'tif');
        saveas(gcf,[savedir1{find(unique(startlocation)==loc)+savedirmodifier},'\','neuron comparingAmplitude_binsize',num2str(binsize),'_',temp,'start',num2str(loc),'.eps'],'epsc');
    end
    
%% single cell plot, after reach obj
    load([fileloc,'\','neuronIndividuals_new.mat']);
    load([fileloc,'\','further_processed_neuron_extraction_final_result.mat']);
    
    temp='C';
    for i=1:length(startlocation)
        neuronIndividuals_new{i}.C=neuronIndividuals_new{i}.C(:,reach_obj_point_all{i}/2:min(floor(away_from_obj_point_all{i})/2-1,length(neuronIndividuals_new{i}.C)));
        neuronIndividuals_new{i}.S=neuronIndividuals_new{i}.S(:,reach_obj_point_all{i}/2:min(floor(away_from_obj_point_all{i})/2-1,length(neuronIndividuals_new{i}.S)));
        neuronIndividuals_new{i}.trace=neuronIndividuals_new{i}.trace(:,reach_obj_point_all{i}/2:min(floor(away_from_obj_point_all{i})/2-1,length(neuronIndividuals_new{i}.trace)));
        neuronIndividuals_new{i}.time=neuronIndividuals_new{i}.time(reach_obj_point_all{i}:min(floor(away_from_obj_point_all{i})-1,length(neuronIndividuals_new{i}.time)));
    end
    
    maxS = max(neuron.C,[],2);
    thresh=3*std(neuron.S,[],2);
    binsize=5;

    for loc=unique(startlocation)
        firingrateAll={};
        countAll={};
        countTimeAll={};
        amplitudeAll={};
        object1=object;
        object1(2)=ROIlist{1}(4)-(object1(2));

        counttt=1;
        for i=1:length(startlocation)
            if startlocation(i)==loc
                [firingrate,count,~,countTime] = calculatingCellSpatialForSingleData_Suoqin(neuronIndividuals_new{i},posAfterObj{i},timeAfterObj{i},ROIlist{i},binsize,1:size(neuronIndividuals_new{i}.C,1),thresh,temp,[],[],[0 1000000]);
                [~,~,~,~,amplitude,amplitude_normalized] = calculatingCellSpatialForSingleData_Suoqin(neuronIndividuals_new{i},posAfterObj{i},timeAfterObj{i},ROIlist{i},binsize,1:size(neuronIndividuals_new{i}.C,1),thresh,temp,[],[],[0 1000000]);  %%%bin size suggests to be 15
%                 [firingrate,count,~,countTime] = calculatingCellSpatialForSingleData_Suoqin(neuronIndividuals_new{i},originalPos{i},originalTime{i},ROIlist{i},binsize,1:size(neuronIndividuals_new{i}.C,1),thresh,temp,[],[],[0 1000000]);
%                 [~,~,~,~,amplitude,amplitude_normalized] = calculatingCellSpatialForSingleData_Suoqin(neuronIndividuals_new{i},originalPos{i},originalTime{i},ROIlist{i},binsize,1:size(neuronIndividuals_new{i}.C,1),thresh,temp,[],[],[0 1000000]);  %%%bin size suggests to be 15
%                 
                firingrateAll(1:size(firingrate,2),counttt)=firingrate;
                countAll(1:size(firingrate,2),counttt)=count;
                for ipp=1:length(firingrate)
                    countTimeAll{ipp,counttt}=countTime;
                end
                amplitudeAll(1:size(firingrate,2),counttt)=amplitude;
                counttt=counttt+1;
            end
        end
        
        for i=1:length(countTimeAll)
            size1=size(countTimeAll{i},1);
            size2=size(countTimeAll{i},2);
        end
        size1m=max(size1);
        size2m=max(size2);
        countTime=zeros(size1m,size2m);
        for i=1:length(countTimeAll)
            for j=1:size(firingrateAll,2)
                countTime=countTime+imresize(countTimeAll{i,j},[size1m,size2m]);
                if ~isempty(firingrateAll{i,j})
                    firingrateAll{i,j}=imresize(firingrateAll{i,j},[size1m,size2m]);
                else
                    firingrateAll{i,j}=zeros(size1m,size2m);
                end
                if ~isempty(countAll{i,j})
                    countAll{i,j}=imresize(countAll{i,j},[size1m,size2m]);
                else
                    countAll{i,j}=zeros(size1m,size2m);
                end
                if ~isempty(amplitudeAll{i,j})
                    amplitudeAll{i,j}=imresize(amplitudeAll{i,j},[size1m,size2m]);
                else
                    amplitudeAll{i,j}=zeros(size1m,size2m);
                end
            end
        end     
        
        firingrateAll_individual=firingrateAll(:,1);
        countAll_individual=countAll(:,1);
        amplitudeAll_individual=amplitudeAll(:,1);
        for j=2:size(firingrateAll,2)
            for i=1:length(countTimeAll)
                if ~isempty(firingrateAll{i,j})
                    firingrateAll_individual{i}=firingrateAll_individual{i}+firingrateAll{i,j};
                end
                if ~isempty(countAll{i,j})
                    countAll_individual{i}=countAll_individual{i}+countAll{i,j};
                end
                if ~isempty(amplitudeAll{i,j})
                    amplitudeAll_individual{i}=amplitudeAll_individual{i}+amplitudeAll{i,j};
                end
            end
        end
        
        plottingFiringBehaviorSpatialForSingleData_adapted_dryland(neuronIndividuals_new(find(startlocation==loc)),posAfterObj(find(startlocation==loc)),timeAfterObj(find(startlocation==loc)),ROIlist{1},countAll_individual,[1:size(neuron.C,1)],thresh,[],[],savedir1(1+savedirmodifier:savedirmodifier+savedirmodifier),find(unique(startlocation)==loc),5,'C',object1,firingrateAll_individual,'firingrate',countTime,loc)
   end

%% vector quantification
tSinceRelease={};
tSinceReleaseRatio={};%ratio to the time reach obj
tPriorObj={};
disSinceRelease={};
disSinceReleaseRatio={};%ratio to the dis reach obj
dis2Obj={};
directDis2Obj={};
%orientation=[]; % to be added when blue light is aviliable

for i=1:length(behavname)
    tSinceRelease{i}=originalTime{i};
    tSinceReleaseRatio{i}=originalTime{i}/originalTime{i}(reach_obj_point_all{i});    
    tPriorObj{i}=originalTime{i}(reach_obj_point_all{i})-originalTime{i};
    
    distancee=[0];
    distance2objLine=[norm(originalPos{i}(reach_obj_point_all{i},:)-originalPos{i}(1,:))];
    for j=2:length(originalPos{i})
        distancee(j)=norm(originalPos{i}(j,:)-originalPos{i}(j-1,:))+distancee(j-1);
        distance2objLine(j)=norm(originalPos{i}(reach_obj_point_all{i},:)-originalPos{i}(j,:));
    end
%     distancee=distancee(2:end);
    
    disSinceRelease{i}=distancee;
    disSinceReleaseRatio{i}=distancee/distancee(reach_obj_point_all{i});
    dis2Obj{i}=abs(distancee(reach_obj_point_all{i})-distancee);
    directDis2Obj{i}=distance2objLine;
end

load([fileloc,'\','neuronIndividuals_new.mat']);
load([fileloc,'\','further_processed_neuron_extraction_final_result.mat']);

temp='C';

% maxS = max(neuron.C,[],2);
% thresh = 0.1*maxS;
thresh=3*std(neuron.S,[],2);
Fs=30;
tau=1/Fs;
nC={};
nS={};
for i=1:length(neuronIndividuals_new)
    downsample=length(neuronIndividuals_new{i}.time)/size(neuronIndividuals_new{i}.C,2);
    nC{i} = interp1(neuronIndividuals_new{i}.time(1:downsample:end),neuronIndividuals_new{i}.C',originalTime{i})'; %%
    nS{i}= interp1(neuronIndividuals_new{i}.time(1:downsample:end),neuronIndividuals_new{i}.S',originalTime{i})'; %%
end

timeRate={};
timeRatioRate={};
timePriorRate={};
disRate={};
disRatioRate={};
dis2ObjRate={};
directDis2ObjRate={};
binNum=50;
% binLength=120;
for i=1:length(neuronIndividuals_new)

    hh1=histogram(tSinceRelease{i},binNum);
    binEdge1=hh1.BinEdges;
    hh2=histogram(tSinceReleaseRatio{i},binNum);
    binEdge2=hh2.BinEdges;
    hh3=histogram(tPriorObj{i},binNum);
    binEdge3=hh3.BinEdges;
    hh4=histogram(disSinceRelease{i},binNum);
    binEdge4=hh4.BinEdges;
    hh5=histogram(disSinceReleaseRatio{i},binNum);
    binEdge5=hh5.BinEdges;
    hh6=histogram(dis2Obj{i},binNum);
    binEdge6=hh6.BinEdges;
    hh7=histogram(directDis2Obj{i},binNum);
    binEdge7=hh7.BinEdges;
%     hh1=sort(tSinceRelease{i});
%     for tk=1:length(hh1)-binLength
%         binEdge1(tk,1)=tk;
%         binEdge1(tk,2)=tk+binLength;
%     end
%     hh2=sort(tSinceReleaseRatio{i});
%     for tk=1:length(hh2)-binLength
%         binEdge2(tk,1)=tk;
%         binEdge2(tk,2)=tk+binLength;
%     end
%     hh3=sort(tPriorObj{i});
%     for tk=1:length(hh3)-binLength
%         binEdge3(tk,1)=tk;
%         binEdge3(tk,2)=tk+binLength;
%     end
%     hh4=sort(disSinceRelease{i});
%     for tk=1:length(hh4)-binLength
%         binEdge4(tk,1)=tk;
%         binEdge4(tk,2)=tk+binLength;
%     end
%     hh5=sort(disSinceReleaseRatio{i});
%     for tk=1:length(hh5)-binLength
%         binEdge5(tk,1)=tk;
%         binEdge5(tk,2)=tk+binLength;
%     end
%     hh6=sort(dis2Obj{i});
%     for tk=1:length(hh6)-binLength
%         binEdge6(tk,1)=tk;
%         binEdge6(tk,2)=tk+binLength;
%     end
%     hh7=sort(directDis2Obj{i});
%     for tk=1:length(hh7)-binLength
%         binEdge7(tk,1)=tk;
%         binEdge7(tk,2)=tk+binLength;
%     end
    close;
%     binNum=size(binEdge1,1);
    
    for j=1:size(neuron.C,1)
        nCt = nC{i}(j,:)>thresh(j); %%
        nSt = nS{i}(j,:)>thresh(j); %%    
%         nCt(nCt<thresh(j))=0; %%
%         nSt(nSt<thresh(j))=0; %%  
        
        for idor3=1:binNum    
            idxx1=double(tSinceRelease{i}>binEdge1(idor3)).*double(tSinceRelease{i}<=binEdge1(idor3+1));
            idxx2=double(tSinceReleaseRatio{i}>binEdge2(idor3)).*double(tSinceReleaseRatio{i}<=binEdge2(idor3+1));
            idxx3=double(tPriorObj{i}>binEdge3(idor3)).*double(tPriorObj{i}<=binEdge3(idor3+1));
            idxx4=double(disSinceRelease{i}>binEdge4(idor3)).*double(disSinceRelease{i}<=binEdge4(idor3+1));
            idxx5=double(disSinceReleaseRatio{i}>binEdge5(idor3)).*double(disSinceReleaseRatio{i}<=binEdge5(idor3+1));
            idxx6=double(dis2Obj{i}>binEdge6(idor3)).*double(dis2Obj{i}<=binEdge6(idor3+1));
            idxx7=double(directDis2Obj{i}>binEdge7(idor3)).*double(directDis2Obj{i}<=binEdge7(idor3+1));            
%             idxx1=double(tSinceRelease{i}>binEdge1(idor3,1)).*double(tSinceRelease{i}<=binEdge1(idor3,2));
%             idxx2=double(tSinceReleaseRatio{i}>binEdge2(idor3,1)).*double(tSinceReleaseRatio{i}<=binEdge2(idor3,2));
%             idxx3=double(tPriorObj{i}>binEdge3(idor3,1)).*double(tPriorObj{i}<=binEdge3(idor3,2));
%             idxx4=double(disSinceRelease{i}>binEdge4(idor3,1)).*double(disSinceRelease{i}<=binEdge4(idor3,2));
%             idxx5=double(disSinceReleaseRatio{i}>binEdge5(idor3,1)).*double(disSinceReleaseRatio{i}<=binEdge5(idor3,2));
%             idxx6=double(dis2Obj{i}>binEdge6(idor3,1)).*double(dis2Obj{i}<=binEdge6(idor3,2));
%             idxx7=double(directDis2Obj{i}>binEdge7(idor3,1)).*double(directDis2Obj{i}<=binEdge7(idor3,2));            

            idxx1=logical(idxx1);
            idxx2=logical(idxx2);
            idxx3=logical(idxx3);
            idxx4=logical(idxx4);
            idxx5=logical(idxx5);
            idxx6=logical(idxx6);
            idxx7=logical(idxx7);
            
            timeRate{i}(idor3,1)=mean(tSinceRelease{i}(idxx1));
            timeRate{i}(idor3,2+j-1)=sum(nSt(idxx1))/(sum(idxx1>0)*tau);
            timeRatioRate{i}(idor3,1)=mean(tSinceReleaseRatio{i}(idxx2));
            timeRatioRate{i}(idor3,2+j-1)=sum(nSt(idxx2))/(sum(idxx2>0)*tau);
            timePriorRate{i}(idor3,1)=mean(tPriorObj{i}(idxx3));
            timePriorRate{i}(idor3,2+j-1)=sum(nSt(idxx3))/(sum(idxx3>0)*tau);    
            disRate{i}(idor3,1)=mean(disSinceRelease{i}(idxx4));
            disRate{i}(idor3,2+j-1)=sum(nSt(idxx4))/(sum(idxx4>0)*tau);   
            disRatioRate{i}(idor3,1)=mean(disSinceReleaseRatio{i}(idxx5));
            disRatioRate{i}(idor3,2+j-1)=sum(nSt(idxx5))/(sum(idxx5>0)*tau);          
            dis2ObjRate{i}(idor3,1)=mean(dis2Obj{i}(idxx6));
            dis2ObjRate{i}(idor3,2+j-1)=sum(nSt(idxx6))/(sum(idxx6>0)*tau);  
            directDis2ObjRate{i}(idor3,1)=mean(directDis2Obj{i}(idxx7));
            directDis2ObjRate{i}(idor3,2+j-1)=sum(nSt(idxx7))/(sum(idxx7>0)*tau);             
        end
    end
end

%plot start
for i=1:length(unique(startlocation))
    for j=1:7
        ax{i,j}=figure(j+(i-1)*7);
    end
end

for loc=unique(startlocation)
    locNum=sum(startlocation==loc);
    ct=1;
    for i=1:length(startlocation)
        if startlocation(i)==loc
            figure(ax{find(unique(startlocation)==loc),1});
            subplot(locNum,1,ct);
            plot(timeRate{i}(:,1),mean(timeRate{i}(:,2:end),2));
            hold on;
            figure(ax{find(unique(startlocation)==loc),2});
            subplot(locNum,1,ct);
            plot(timeRatioRate{i}(:,1),mean(timeRatioRate{i}(:,2:end),2));
            hold on;
            figure(ax{find(unique(startlocation)==loc),3});
            subplot(locNum,1,ct);
            plot(timePriorRate{i}(timePriorRate{i}(:,1)>0,1),mean(timePriorRate{i}(timePriorRate{i}(:,1)>0,2:end),2));
            hold on;     
            figure(ax{find(unique(startlocation)==loc),4});
            subplot(locNum,1,ct);
            plot(disRate{i}(:,1),mean(disRate{i}(:,2:end),2));
            hold on; 
            figure(ax{find(unique(startlocation)==loc),5});
            subplot(locNum,1,ct);
            plot(disRatioRate{i}(:,1),mean(disRatioRate{i}(:,2:end),2));
            hold on;           
            figure(ax{find(unique(startlocation)==loc),6});
            subplot(locNum,1,ct);
            plot(dis2ObjRate{i}(:,1),mean(dis2ObjRate{i}(:,2:end),2));
            hold on;       
            figure(ax{find(unique(startlocation)==loc),7});
            subplot(locNum,1,ct);
            plot(directDis2ObjRate{i}(:,1),mean(directDis2ObjRate{i}(:,2:end),2));
            hold on; 
            ct=ct+1;
        end
    end
end

activityVecName={'time since release','time since release (Ratio)','time prior reach Obj','distance since release','distance since release (Ratio)','travel distance to Obj','straight distance to Obj'}; 
mkdir([savedir,'\','Activity Rate Vector']);
for i=1:length(unique(startlocation))
    ust=unique(startlocation);
    locNum=sum(startlocation==loc);
    for j=1:7
        figure(ax{i,j});
        title(['start ',num2str(ust(i)),' ',activityVecName{j}]);
        xlabel(activityVecName{j});
        ylabel('event rate');
        saveas(gcf,[savedir,'\','Activity Rate Vector','\','start',num2str(ust(i)),activityVecName{j},'.fig'],'fig');
        saveas(gcf,[savedir,'\','Activity Rate Vector','\','start',num2str(ust(i)),activityVecName{j},'.eps'],'epsc');
        saveas(gcf,[savedir,'\','Activity Rate Vector','\','start',num2str(ust(i)),activityVecName{j},'.tif'],'tif');        
    end
end
close all
