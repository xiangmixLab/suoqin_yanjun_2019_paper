%vname gen
for i=1:length(orilocation)
vname{i}=orilocation{i};
labelll=findstr(vname{i},'\');
vname{i}(1:labelll(3))=[];
vname{i}(findstr(vname{i},'\'))='_';
end

for i=1:length(orilocation)
    cd(orilocation{i})
    v=VideoReader('behavCam1_1.avi');
    count=1;
    vid=[];
    behavROI=[];
    while hasFrame(v)
        t=readFrame(v);
        vid(:,:,count)=squeeze(t(:,:,3));
        count=count+1;
    end
   
    medianvid=uint8(median(vid,3));
    medianvidlist{i}=medianvid;
end
objlist={};
%ROI
for i=1:length(orilocation)
     medianvidt=medianvidlist{i};
     medianvid=medianvidlist{i};
%     medianvid_reshape=reshape(medianvid,size(medianvid,1)*size(medianvid,2),1);
%     group=kmeans(medianvid_reshape,3);
%     group_reshape=reshape(group,size(medianvid,1),size(medianvid,2));
%     unique_group_area=[sum(group==1),sum(group==2),sum(group==3)];
%     behavROI=group_reshape==find(unique_group_area==max(unique_group_area));
%     behavROI=imfill(behavROI,'holes');
    medianvid=imadjust(medianvid,[0.07,0.08]);
    behavROI=medianvid<255;
    
    se=strel('disk',15);
    se1=strel('disk',5);
    behavROI=imopen(behavROI,se);
    behavROI=imfill(behavROI,'holes');
    behavROI=bwareaopen(behavROI,10000);
    behavROI=imclose(behavROI,se1); 
    
    behavROI=imerode(behavROI,se1);
    behavROI=bwareaopen(behavROI,10000);
    se2=strel('disk',30);
    behavROI=imfill(behavROI,'holes');
    behavROI=imclose(behavROI,se2);
    medianvidt=medianvidt.*uint8(behavROI);
%     medianvid=imclose(medianvid,se);
    imagesc(medianvidt);
    ROIlist{i}=logical(behavROI);
end
%obj
for i=1:length(orilocation)
     medianvidt=medianvidlist{i};
     medianvid=medianvidlist{i};
    medianvid=imadjust(medianvid,[0.07,0.08]);
    medianvidt=medianvidt.*uint8(behavROI);
    imagesc(medianvidt);
    objloct=[];
    [objloct(:,1),objloct(:,2)]=getpts;   
%     objloc=objloct(1,:);
    objloc=objloct;
    if behavROI(round(objloc(2)),round(objloc(1)))==0
        objloc=[0 0];
    end
    
    ROIlist{i}=logical(behavROI);
    objlist{i}=objloc;
    disp(['finish',num2str(i)]);
end
    
    
    
%     [~,~,bw_all,behav]=Miniscope_behav_extraction_auto_dry_land(orilocation,vname,'behavCam1_',[1:length(orilocation)],'ROIlist',ROIlist,'List',objlist);