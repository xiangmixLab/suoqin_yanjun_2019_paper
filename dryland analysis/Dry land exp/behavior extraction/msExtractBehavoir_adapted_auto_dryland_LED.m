function behav = msExtractBehavoir_adapted_auto_dryland_LED(behav, trackLength,startframe)
%MSEXTRACTBEHAVOIR Summary of this function goes here
%   Detailed explanation goes here
    ROI = uint16([floor(behav.ROI(1)) floor(behav.ROI(1))+behav.ROI(3) floor(behav.ROI(2)) floor(behav.ROI(2))+behav.ROI(4)]);
%% Removes background
% takes the median of each pixel from 100 frames spread out across the video
    numFramesUsed = 100;
    if behav.numFrames >  numFramesUsed
        backgroundFrames = ceil(linspace(1,behav.numFrames,numFramesUsed));
    else
        backgroundFrames = linspace(1,behav.numFrames,behav.numFrames);
    end
    
    frame = uint8(zeros(behav.height,behav.width,3,length(backgroundFrames)));
    for index=1:length(backgroundFrames)
        if (mod(index,10)==0)
            display(['Reading in video for background subtraction. ' num2str(index/length(backgroundFrames)*100) '% done'])
        end
        framet=uint8(msReadFrame(behav,backgroundFrames(index),false,false,false));
        framet=imresize(framet,[behav.height,behav.width]);
        frame(:,:,:,index) = framet;
    end
    frame = frame(ROI(3):ROI(4),ROI(1):ROI(2),:,:);
    background = median(frame,4);
    figure(2)
    imshow(background,'InitialMagnification','fit')
    title('Background. Make sure the mouse does not show up')
%%
    position = nan(behav.numFrames, 2);
    positionblue = nan(behav.numFrames,2);
    for frameNum=1:behav.numFrames
        %% using centroid to replace rgb/hsv value added 022119
        position(frameNum,[1 2])=[behav.rcentroid(frameNum,1)-behav.ROI(1),behav.rcentroid(frameNum,2)-behav.ROI(2)];
        positionblue(frameNum,[1 2])=[behav.bcentroid(frameNum,1)-behav.ROI(1),behav.bcentroid(frameNum,2)-behav.ROI(2)];
    end
    position(1:startframe,:)=repmat([nan nan],length(1:startframe),1);% first frame might be the recording from last run
    position(end,:)=[nan nan];% last frame might also be weird
    position = position*trackLength/behav.ROI(3);
%     idx_nonNan=find(~isnan(position(:,1)));
%     position(1,:)=position(idx_nonNan(2),:);
    time = behav.time(~isnan(position(:,1)));
    position = position(~isnan(position(:,1)),:);
    
%     %  add new code to debug the interpolation problem by Suoqin  
%     temp = find(diff(time)==0);
%     position(temp,:) = [];
%     time(temp) = [];
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% by Suoqin   

    position = interp1(time,position,behav.time,'linear');
    
%     position(:,2)=263-position(:,2);
    dt = median(diff(behav.time/1000)); 
%     position = smoothts(position','b',ceil(1/dt))';
%    [position(:,1),position(:,2)]=smooth_contours(position(:,1), position(:,2), ceil(1/dt));
    position = smoothdata(position,'SmoothingFactor',0.05);
    behav.position = position;
    
    dx = [0; diff(position(:,1))];
    dy = [0; diff(position(:,2))];
    
    dl = sqrt((dx).^2+(dy).^2);
    behav.distancered = sum(dl);
    behav.speed = sqrt((dx).^2+(dy).^2)/dt;
%     behav.speed = smoothts(behav.speed','b',ceil(1/dt));
    behav.dt = dt;
    behav.trackLength = trackLength;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% blue LED
    positionblue(1:startframe,:)=repmat([nan nan],length(1:startframe),1);% first 5 frame might be the recording from last run
    positionblue(end,:)=[nan nan];% last 5 frame might also be weird
    positionblue = positionblue*trackLength/behav.ROI(3);
%     idx_nonNan=find(~isnan(positionblue(:,1)));
%     positionblue(1,:)=positionblue(idx_nonNan(2),:);
    time = behav.time(~isnan(positionblue(:,1)));
    positionblue = positionblue(~isnan(positionblue(:,1)),:);  

    positionblue = interp1(time,positionblue,behav.time,'linear');
%     positionblue(:,2)=263-positionblue(:,2);
    dt = median(diff(behav.time/1000)); 
%     positionblue = smoothts(positionblue','b',ceil(1/dt))';
%    [positionblue(:,1) positionblue(:,2)]=smooth_contours(positionblue(:,1), positionblue(:,2), ceil(1/dt));
    positionblue = smoothdata(positionblue,'SmoothingFactor',0.05);   
    behav.positionblue = positionblue;
    
    dx = [0; diff(positionblue(:,1))];
    dy = [0; diff(positionblue(:,2))];
    
    dl = sqrt((dx).^2+(dy).^2);
    behav.distanceblue = sum(dl);
    behav.speedblue = sqrt((dx).^2+(dy).^2)/dt;
%     behav.speedblue = smoothts(behav.speedblue','b',ceil(1/dt));
    behav.dt = dt;
    behav.trackLength = trackLength;    
end


