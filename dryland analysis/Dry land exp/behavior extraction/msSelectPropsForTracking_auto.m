
function [behav,frame,frame1_bg] = msSelectPropsForTracking_auto(behav,ROImethod,ListedROI,o1)
%MSSELECTPROPSFORTRACKING Summary of this function goes here
%   Detailed explanation goes here
    userInput = 'N';
    green = cat(3, zeros(behav.height,behav.width), ...
    ones(behav.height,behav.width), ...
    zeros(behav.height,behav.width));

    area_ratio_thres=0.0007;
    frame = double(msReadFrame(behav,round(behav.numFrames/2),false,false,false))/255;
    
    tic
    count=1;
    frame1=zeros(size(frame,1),size(frame,2),3,round(behav.numFrames/(10)));
    for i=1:10:behav.numFrames
        t=double(msReadFrame(behav,i,false,false,false))/255;
        if size(t,1)~=size(frame,1)||size(t,2)~=size(frame,2)
            t=imresize(t,[size(frame,1),size(frame,2)]);
        end
        frame1(:,:,:,count)=t;
        count=count+1;
    end
    frame1_bg=median(frame1,4);
    
    toc;
    clear frame1;
    %% ROI generation
    if isequal(ROImethod,'autoROI')
        frame1_bg_g=rgb2gray(frame1_bg);
        Edge_frame1_bg_g_1 = edge(frame1_bg_g,'canny');
        se1 = strel('disk',10);
        Edge_frame1_bg_g_1=imclose(Edge_frame1_bg_g_1,se1);
        frame1_bg_g=adapthisteq(frame1_bg_g);
        se = strel('disk',40);
        Ie = imerode(frame1_bg_g,se);
        frame1_bg_g_1 = imreconstruct(Ie,frame1_bg_g);
        frame1_bg_resize=reshape(frame1_bg_g_1,size(frame1_bg,1)*size(frame1_bg,2),1);
        kindex=kmeans(frame1_bg_resize,3);
        kindex_resize=reshape(kindex,size(frame1_bg,1),size(frame1_bg,2));

        stats=regionprops(kindex_resize,'Area'); 
        area_all=[stats.Area];
        kindex_resize(kindex_resize~=find(area_all==max(area_all)))=0;
        kindex_resize(Edge_frame1_bg_g_1==1)=0;
        BWao = bwareaopen(kindex_resize,6000);
        se2 = strel('disk',20);
        BWao1 = imclose(BWao,se2);
        stats1=regionprops(BWao1,'Area', 'BoundingBox'); 
        rect=stats1.BoundingBox;
        behav.ROI = rect;
    else if isequal(ROImethod,'manualROI')        
            while (strcmp(userInput,'N'))
                clf
                imshow(frame,'InitialMagnification','fit');
                hold on
                display('Select ROI');
                rect = getrect(); 

                behav.ROI = rect; %uint16([rect(1) rect(1)+rect(3) rect(2) rect(2)+rect(4)]);
                rectangle('Position',rect,'LineWidth',2);
                hold off
                userInput = upper(input('Keep ROI? (Y/N)','s'));
            end
    else if isequal(ROImethod,'ROIlist')&&~isempty(ListedROI)
                rect=ListedROI{o1};
                if size(rect,1)==1&&size(rect,2)==4
                    behav.ROI = ListedROI{o1};
                end
                if size(rect,1)==size(frame1_bg,1)&&size(rect,2)==size(frame1_bg,2)
                    statss=regionprops(rect,'BoundingBox');
                    behav.ROI=statss.BoundingBox;
                end
         end
        end           
    end
    
    %% find global maximum
    tic;
    max_framer=[];
    max_frameb=[];
    count=1;
    for i=1:10:behav.numFrames
        frame1t=double(msReadFrame(behav,i,false,false,false))/255;
        if size(frame1t,1)~=size(frame,1)||size(frame1t,2)~=size(frame,2)
            frame1t=imresize(frame1t,[size(frame,1),size(frame,2)]);
        end
        t=squeeze(frame1t)-frame1_bg;
        if size(rect,1)==1&&size(rect,2)==4
            t_lab=rgb2lab(t(round(rect(2)):round(rect(2)+rect(4)),round(rect(1)):round(rect(1)+rect(3)),:));
        else
            t_lab=rgb2lab(t.*double(rect));
        end
        max_framer(count)=max(max(t_lab(:,:,2)));
        max_frameb(count)=max(max(t(:,:,3)));
        count=count+1;
    end
    max_framer1=median(max_framer);
    max_frameb1=median(max_frameb);
    toc;
    
    %% red and blue LED
    se=strel('disk',3);
    se1=strel('disk',5);
    tic;
    
%     rcentroid=nan(behav.numFrames,2);
%     bcentroid=nan(behav.numFrames,2);
%     hsvLevelsr=nan(behav.numFrames,6);
%     rgbLevelsr=nan(behav.numFrames,6);
%     hsvLevelsb=nan(behav.numFrames,6);
%     rgbLevelsb=nan(behav.numFrames,6);
    parfor i=2:behav.numFrames
        [hsvLevelsr(i,:),rgbLevelsr(i,:),hsvLevelsb(i,:),rgbLevelsb(i,:),rcentroid(i,:),bcentroid(i,:)]=LED_position_extract(behav,i,frame1_bg,rect,se,se1,max_framer1,max_frameb1);       
        toc;        
    end
    
    [y,x]=find(isnan(hsvLevelsr));
    hsvLevelsr(y,:)=repmat([-1 -1 -1 -1 -1 -1],length(y),1);
    behav.hsvLevelred = hsvLevelsr;
    [y,x]=find(isnan(rgbLevelsr));
    rgbLevelsr(y,:)=repmat([-1 -1 -1 -1 -1 -1],length(y),1);
    behav.rgbLevelred = rgbLevelsr;
    [y,x]=find(isnan(rcentroid));
    rcentroid(y,:)=repmat([nan nan],length(y),1);
    behav.rcentroid=rcentroid;
    
    [y,x]=find(isnan(hsvLevelsb));
    hsvLevelsb(y,:)=repmat([-1 -1 -1 -1 -1 -1],length(y),1);
    behav.hsvLevelblue = hsvLevelsb;
    [y,x]=find(isnan(rgbLevelsb));
    rgbLevelsb(y,:)=repmat([-1 -1 -1 -1 -1 -1],length(y),1);
    behav.rgbLevelblue = rgbLevelsb;
    [y,x]=find(isnan(bcentroid));
    bcentroid(y,:)=repmat([nan nan],length(y),1);
    behav.bcentroid=bcentroid;
    

