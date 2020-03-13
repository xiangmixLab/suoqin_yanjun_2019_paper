%ROI:[x,y,lx,ly]
%object coordination: left bottom [0,0]
%ROImethod: autoROI, manualROI, ROIlist
%objInput: Y,N,List
function [behavfile_list,camera_id_list]=Miniscope_behav_extraction_auto_dryland_LED(orilocation,vname,behavVidprefix,mrange,ROImethod,ROIlist,objInput,objlist,trackLength,startframe)
tic;
behavfile_list={};
camera_id_list=[];

for ol=mrange
    
behavloc=orilocation{ol};
if ~isempty(behavloc)
cd(behavloc);
%% Generate behav.m
behav = msGenerateVideoObj_auto(pwd,behavVidprefix);

%% Select ROI and HSV for tracking
[behav,framest,bg] = msSelectPropsForTracking_auto(behav,ROImethod,ROIlist,ol); 
size(framest)
%% TRACK LENGTH
% trackLength = 210;%cm  % box experiments
% trackLength = 280;%cm  % box experiments
% trackLength = 536;%cm  % linear track experiments %% used scaled track length
% trackLength = 420;%cm  % FN box experiments


%% Select Objects
% display('Select Object');
% userInput = upper(input('Select Object? (Y/N)','s'));
if (strcmp(objInput,'Y'))
    behav = SelectObject(behav,trackLength);% as we are using ginput in image, the object location is in xy axis ((0,0) left bottom)
else
    behav.object=[0,0];
end

if (strcmp(objInput,'auto'))
    see=strel('disk',15);
    if length(size(bg))>2
        bg1=rgb2gray(bg);
    else
        bg1=bg;
    end
    bg1(:,[1:10,end-10:end])=0;
    bg1([1:10,end-10:end],:)=0;
    bg1=imclose(bg1,see);
    Centroid_obj=[0,0];
    for thres_obj=0:0.01:1
        bgt=bg1>thres_obj;
        bgt=bwareaopen(bgt,50);
        statss=regionprops(logical(bgt));
        Area_all=[statss.Area];
        if length(Area_all)<5&&sum(Area_all)/(size(bg1,1)*size(bg1,2))<0.02
            thresh_obj=thres_obj;
            for cen_num=1:length(Area_all)
                Centroid_obj(cen_num,:)=statss(cen_num).Centroid;
            end
            break;
        end
    end
    behav.object=Centroid_obj;
end
co=[];
if (strcmp(objInput,'List'))&&~isempty(objlist)% object: just maintain the coordination selected from imagej or ginput
    if sum(objlist{ol})>0
    for lkk=1:size(objlist{ol},1)
        co(lkk,:)=[objlist{ol}(lkk,1)-behav.ROI(1),behav.ROI(2)+behav.ROI(4)-objlist{ol}(lkk,2)];
    end
    else
        co=[0,0];
    end
    behav.object=co;
end

%%% new code from Tristan 10/06/2016%%%%
%% Extract position
% important variable, which is set according to the video trace
 
%correct time variable  %%05/2017 modification 
d=diff(behav.time);
t=find(d<=0);
while ~isempty(t)
behav.time(t+1)=behav.time(t)+1;
d=diff(behav.time);
t=find(d<=0);
end
behav = msExtractBehavoir_adapted_auto_dryland_LED(behav, trackLength,startframe(ol)); 
% 
% if ~isempty(isnan(behav.position))
%     nanind1=find(isnan(behav.position(:,1)));
%     behav.position(nanind1,1)=behav.position(max(nanind1)+1,1);
%     nanind2=find(isnan(behav.position(:,2)));
%     behav.position(nanind2,2)=behav.position(max(nanind2)+1,2);
% end
% if ~isempty(isnan(behav.positionblue))
%     nanind1=find(isnan(behav.positionblue(:,1)));
%     behav.positionblue(nanind1,1)=behav.positionblue(max(nanind1)+1,1);
%     nanind2=find(isnan(behav.positionblue(:,2)));
%     behav.positionblue(nanind2,2)=behav.positionblue(max(nanind2)+1,2);
% end

% ROI3=behav.ROI(3);
% behav.ROI=behav.ROI*trackLength/ROI3;
% 
% % adjust object pos to real box size
% behav.object=behav.object*trackLength/ROI3;

%% velocity calculation to see stationary points
velocityx=velocity_kevin_cal(behav.time,behav.position(:,1));
velocityy=velocity_kevin_cal(behav.time,behav.position(:,2));
velocity=(velocityx.^2+velocityy.^2).^0.5;
velosmallthan3=velocity<(3*behav.ROI(1,3)/behav.trackLength);
% behav.time(velosmallthan3)=[];
% behav.position(velosmallthan3)=[];

%% mouse head orientation
% mouseheading=[behav.position(:,1)-behav.positionblue(:,1),behav.position(:,2)-behav.positionblue(:,2)];
% mouse_object_dirction={};
% for i=1:size(behav.object,1)
%     mouse_object_dirction{i}=[repmat(behav.object(i,1),length(behav.positionblue(:,1)),1)-behav.positionblue(:,1),repmat(behav.object(i,2),length(behav.positionblue(:,2)),1)-behav.positionblue(:,2)];
% end
% if_mouse_head_toward_object=zeros(size(mouseheading,1),length(mouse_object_dirction));
% mouse_head_toward_object_angle=zeros(size(mouseheading,1),length(mouse_object_dirction));
% for i=1:size(behav.object,1)
%     for j=1:size(mouseheading,1)
%         ThetaInDegrees = atan2d(norm(cross([mouseheading(j,:) 0],[mouse_object_dirction{i}(j,:) 0])),dot([mouseheading(j,:) 0],[mouse_object_dirction{i}(j,:) 0]));
%         mouse_head_toward_object_angle(j,i)=ThetaInDegrees;
%     end 
% end
% for i=1:size(behav.object,1)
%     for j=1:size(mouseheading,1)-15
%         if abs(mean(mouse_head_toward_object_angle(j:j+14,i)))<=30
%             if_mouse_head_toward_object(j,i)=1;
%         end
%     end 
% end
% if size(if_mouse_head_toward_object,2)>1
% for j=1:size(mouseheading,1)
%         if if_mouse_head_toward_object(j,1)==1&&if_mouse_head_toward_object(j,2)==1
%             L1=((behav.object(1,1)-behav.position(j,1))^2+(behav.object(1,2)-behav.position(j,2))^2)^0.5;
%             L2=((behav.object(2,1)-behav.position(j,2))^2+(behav.object(2,2)-behav.position(j,2))^2)^0.5;
%             if L1>L2
%                 if_mouse_head_toward_object(j,2)=0;
%             end
%             if L1<L2
%                 if_mouse_head_toward_object(j,1)=0;
%             end                
%         end
% end 
% end
% 
% mouse_head_object_cal.mouseheading=mouseheading;
% mouse_head_object_cal.mouse_object_dirction=mouse_object_dirction;
% mouse_head_object_cal.mouse_head_toward_object_angle=mouse_head_toward_object_angle;
% mouse_head_object_cal.if_mouse_head_toward_object=if_mouse_head_toward_object;

% filename=dir();
% for i=1:10
%     if ~isempty(strfind(filename(i).name,'NormCorre'))
%         filename1=filename(i).name(11:end-4);
%     end
% end

filename1=vname{ol};

%% Check behav trace using the lines below. Include figures into your notes
figure; hist(diff(behav.time))
saveas(gcf,'behavtimeinterval.tif');
saveas(gcf,'behavtimeinterval.fig');
figure; 
imshow(framest);hold on;
rectangle('Position',behav.ROI,'LineWidth',2);
plot(behav.position(:,1)*behav.ROI(3)/trackLength+behav.ROI(1),behav.position(:,2)*behav.ROI(3)/trackLength+behav.ROI(2))
if sum(behav.object)>0
    plot(behav.object(:,1)+behav.ROI(1),behav.ROI(4)+behav.ROI(2)-behav.object(:,2),'.','MarkerSize',20)
end
axis image
saveas(gcf,'behavtrack_auto_red.png');
saveas(gcf,'behavtrack_auto_red.fig');
figure; 
imshow(framest);hold on;
rectangle('Position',behav.ROI,'LineWidth',2);
plot(behav.positionblue(:,1)*behav.ROI(3)/trackLength+behav.ROI(1),behav.positionblue(:,2)*behav.ROI(3)/trackLength+behav.ROI(2))
if sum(behav.object)>0
    plot(behav.object(:,1)+behav.ROI(1),behav.ROI(4)+behav.ROI(2)-behav.object(:,2),'.','MarkerSize',20)
end
axis image
saveas(gcf,'behavtrack_auto_blue.tif');
saveas(gcf,'behavtrack_auto_blue.fig');
close all;

% adjust ROI to real box size
ROI3=behav.ROI(3);
behav.ROI=behav.ROI*behav.trackLength/ROI3;
behav.ROI3=ROI3;
% adjust object pos to real box size
behav.object=behav.object*behav.trackLength/ROI3;

behavfile_list{ol}=[filename1,'_','auto_Behav.mat'];
camera_id_list(ol,:)=[behav.msCam_num,behav.camNumber];
save([filename1,'_','Behav.mat'], 'behav','velosmallthan3')
toc;
disp(['finish ',num2str(ol)])
end
end