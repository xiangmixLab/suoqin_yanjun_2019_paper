
function [behav,frame1_bg] = msExtractBehavior_adapted_auto_dry_land(behav,ROImethod,ListedROI,o1)
%MSSELECTPROPSFORTRACKING Summary of this function goes here
%   Detailed explanation goes here
    userInput = 'N';
    green = cat(3, zeros(behav.height,behav.width), ...
    ones(behav.height,behav.width), ...
    zeros(behav.height,behav.width));


    frame = double(msReadFrame(behav,round(behav.numFrames/2)+100,false,false,false))/255;
    
    tic
    count=1;
    frame1=zeros(size(frame,1),size(frame,2),3,round(behav.numFrames/10));
    bwr=zeros(size(frame,1),size(frame,2));
    for i=1:10:behav.numFrames
        frame1(:,:,:,count)=double(msReadFrame(behav,i,false,false,false))/255;
        bwr=bwr+double(imbinarize(frame1(:,:,1,count),0.5));
        count=count+1;
    end
    bwr=logical(bwr);
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
    %     if isequal(cuthalf,'half')
    %         statbwr=regionprops(bwr);
    %         centroid=statbwr.Centroid;
    %         if centroid(1)<=0.5*size(frame,2)
    %             kindex_resize=kindex_resize(:,1:round(0.5*size(frame,2)));
    %         else
    %             kindex_resize=kindex_resize(:,round(0.5*size(frame,2)):end);
    %         end
    %     end
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
                behav.ROI = ListedROI{o1}; %uint16([rect(1) rect(1)+rect(3) rect(2) rect(2)+rect(4)]);
                rect=ListedROI{o1};
        end
        end           
    end
    %% ear and tail
    
    for i=1:1:behav.numFrames
        tempf=double(msReadFrame(behav,i,false,false,false))/255;
        tempf=tempf-frame1_bg;
        frame2=tempf(round(rect(2)):round(rect(2)+rect(4)),round(rect(1)):round(rect(1)+rect(3)),:);
        frame2_gr=rgb2gray(squeeze((frame2)));

        tic;

        tempf2=zeros(size(tempf,1),size(tempf,2));
        tempf2(round(rect(2)):round(rect(2)+rect(4)),round(rect(1)):round(rect(1)+rect(3)))=frame2_gr;
        bw=imbinarize(tempf2);
        bw1=~bwareaopen(bw,1000);
        bw2=bwareaopen(bw,9);
        bw=bw(bw1.*bw2);
        bwl=bwlabel(bw);
        bwstats=regionprops(bwl,'Area','MajorAxisLength','MinorAxisLength');
        AxisRatio=[bwstats.MajorAxisLength]./[bwstats.MinorAxisLength];
        hsvLevelsb(i,:) = [H+[-.2 .2] S+[-.2 .2] V+[-.2 .2]];
        rgbLevelsb(i,:) = [R+[-.2 .2] G+[-.2 .2] B+[-.2 .2]]; 
    end
    
    [y,x]=find(isnan(hsvLevelsr));
    hsvLevelsr(y,:)=repmat([-1 -1 -1 -1 -1 -1],length(y),1);
    behav.hsvLevelred = hsvLevelsr;
    [y,x]=find(isnan(rgbLevelsr));
    rgbLevelsr(y,:)=repmat([-1 -1 -1 -1 -1 -1],length(y),1);
    behav.rgbLevelred = rgbLevelsr;
    
    [y,x]=find(isnan(hsvLevelsb));
    hsvLevelsb(y,:)=repmat([-1 -1 -1 -1 -1 -1],length(y),1);
    behav.hsvLevelblue = hsvLevelsb;
    [y,x]=find(isnan(rgbLevelsb));
    rgbLevelsb(y,:)=repmat([-1 -1 -1 -1 -1 -1],length(y),1);
    behav.rgbLevelblue = rgbLevelsb;

        position = position*trackLength/behav.ROI(3);
    time = behav.time(~isnan(position(:,1)));
    position = position(~isnan(position(:,1)),:);
    
%     %  add new code to debug the interpolation problem by Suoqin  
%     temp = find(diff(time)==0);
%     position(temp,:) = [];
%     time(temp) = [];
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% by Suoqin   

    position = interp1(time,position,behav.time);
%     position(:,2)=263-position(:,2);
    dt = median(diff(behav.time/1000)); 
    position = smoothts(position','b',ceil(1/dt))';
    behav.position = position;
    
    dx = [0; diff(position(:,1))];
    dy = [0; diff(position(:,2))];
    
    dl = sqrt((dx).^2+(dy).^2);
    behav.distancered = sum(dl);
    behav.speed = sqrt((dx).^2+(dy).^2)/dt;
    behav.speed = smoothts(behav.speed','b',ceil(1/dt));
    behav.dt = dt;
    behav.trackLength = trackLength;
    
    
    positionblue = positionblue*trackLength/behav.ROI(3);
    time = behav.time(~isnan(positionblue(:,1)));
    positionblue = positionblue(~isnan(positionblue(:,1)),:);  

    positionblue = interp1(time,positionblue,behav.time);
%     positionblue(:,2)=263-positionblue(:,2);
    dt = median(diff(behav.time/1000)); 
    positionblue = smoothts(positionblue','b',ceil(1/dt))';
    behav.positionblue = positionblue;
    
    dx = [0; diff(positionblue(:,1))];
    dy = [0; diff(positionblue(:,2))];
    
    dl = sqrt((dx).^2+(dy).^2);
    behav.distanceblue = sum(dl);
    behav.speedblue = sqrt((dx).^2+(dy).^2)/dt;
    behav.speedblue = smoothts(behav.speedblue','b',ceil(1/dt));
    behav.dt = dt;
    behav.trackLength = trackLength;    

