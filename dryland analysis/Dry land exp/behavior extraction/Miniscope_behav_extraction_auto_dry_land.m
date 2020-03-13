%ROI:[x,y,lx,ly]
%object coordination: left bottom [0,0]
%ROImethod: autoROI, manualROI, ROIlist
%objInput: Y,N,List
function [behavfile_list,camera_id_list,bw_all,behav]=Miniscope_behav_extraction_auto_dry_land(orilocation,vname,behavVidprefix,mrange,ROImethod,ROIlist,objInput,objlist)

behavfile_list={};
camera_id_list=[];
for ol=mrange
    
behavloc=orilocation{ol};
cd(behavloc);

%% Generate behav.m
behav = msGenerateVideoObj_auto(pwd,behavVidprefix);

%% TRACK LENGTH
% trackLength = 280;%cm  % box experiments
% trackLength = 536;%cm  % linear track experiments %% used scaled track length
% trackLength = 420;%cm  % FN box experiments
trackLength = behav.vidObj{1,1}.Width;%pixel

%% select objects
if (strcmp(objInput,'Y'))
    behav = SelectObject(behav,trackLength);% as we are using ginput in image, the object location is in xy axis ((0,0) left bottom)
else
    behav.object=[0,0];
end

if (strcmp(objInput,'List'))&&~isempty(objlist)
    behav.object=objlist{ol};
end


%% Extract position
[behav,frame,bw_all] = msExtractBehavior_adapted_auto_dry_land(behav,ROImethod,ROIlist,ol,trackLength); 

%% 

%% velocity calculation to see stationary points
velocityx=velocity_kevin_cal(behav.time,behav.positionEar(:,1));
velocityy=velocity_kevin_cal(behav.time,behav.positionEar(:,2));
velocity=(velocityx.^2+velocityy.^2).^0.5;
velosmallthan3=velocity<(3*behav.ROI(1,3)/behav.trackLength);

filename1=vname{ol};

save([filename1,'_','Behav.mat'], 'behav','velosmallthan3')

%% Check behav trace using the lines below. Include figures into your notes
figure; hist(diff(behav.time))
saveas(gcf,'behavtimeinterval.tif');
saveas(gcf,'behavtimeinterval.fig');
figure; 
imshow(frame);hold on;
rectangle('Position',behav.ROI,'LineWidth',2);
plot(behav.positionEar(:,1),behav.positionEar(:,2),'c-','lineWidth',1.5);
if sum(behav.object)>0
    plot(behav.object(:,1),behav.object(:,2),'r.','MarkerSize',30)
end
axis image
saveas(gcf,'behavtrack_Ear.tif');
saveas(gcf,'behavtrack_Ear.eps','epsc');
saveas(gcf,'behavtrack_Ear.fig');
figure; 
imshow(frame);hold on;
rectangle('Position',behav.ROI,'LineWidth',2);
plot(behav.positionTail(:,1),behav.positionTail(:,2),'c-','lineWidth',1.5)
if sum(behav.object)>0
    plot(behav.object(:,1),behav.object(:,2),'.','MarkerSize',30)
end
axis image
saveas(gcf,'behavtrack_Tail.tif');
saveas(gcf,'behavtrack_Tail.eps','epsc');
saveas(gcf,'behavtrack_Tail.fig');
close all;

behavfile_list{ol}=[filename1,'_','Behav.mat'];
camera_id_list(ol,:)=[behav.msCam_num,behav.camNumber];
disp(['finish ',num2str(ol)]);

%% quick illustration
% for p=229:length(orilocation)
% cd(orilocation{p});
% load([vname{p},'_Behav.mat']);
% v = VideoWriter(['mouseDirIllustration.avi']);
% open(v);
% tempf=[];
% tempf=uint8(tempf);
%  for i=1:1:behav.numFrames
% tempf(:,:,i)=rgb2gray(uint8(msReadFrame(behav,i,false,false,false)));
% end
% for i=1:size(tempf,3)
% imshow(tempf(:,:,i));hold on
% quiver(behav.positionTail(i,1),behav.positionTail(i,2),behav.positionEar(i,1)-behav.positionTail(i,1),behav.positionEar(i,2)-behav.positionTail(i,2),'MaxHeadSize',1000);
% title(num2str(i));
% drawnow;
% frame=getframe(gcf);
% writeVideo(v,frame);
% clf
% end
% close(v);
% disp(['finish ',num2str(p)]);
% end
end