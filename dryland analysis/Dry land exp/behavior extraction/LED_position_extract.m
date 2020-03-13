function [hsvLevelsr,rgbLevelsr,hsvLevelsb,rgbLevelsb,rcentroid,bcentroid]=LED_position_extract(behav,i,frame1_bg,rect,se,se1,max_framer,max_frameb)       
        tempf=double(msReadFrame(behav,i,false,false,false))/255;
        if size(tempf,1)~=size(frame1_bg,1)||size(tempf,2)~=size(frame1_bg,2)
            tempf=imresize(tempf,[size(frame1_bg,1),size(frame1_bg,2)]);
        end
        if length(size(frame1_bg))>3
            tempf=tempf-squeeze(frame1_bg(:,:,:,min(ceil(i/40),floor(behav.numFrames/40))));
        else
            tempf=tempf-frame1_bg;
        end
        if size(rect,1)==1&&size(rect,2)==4
            frame2=tempf(round(rect(2)):round(rect(2)+rect(4)),round(rect(1)):round(rect(1)+rect(3)),:);           
        end
        if size(rect,1)==size(tempf,1)&&size(rect,2)==size(tempf,2)
            frame2t=tempf.*double(imdilate(rect,se1));
            frame2=frame2t(floor(behav.ROI(2)):floor(behav.ROI(2)+behav.ROI(4)),floor(behav.ROI(1)):floor(behav.ROI(1)+behav.ROI(3)),:);
        end
        
        frame2_r=squeeze(frame2(:,:,1));
        frame2_b=squeeze(frame2(:,:,3));

        r_b_diff=abs(sum(sum(frame2_r-frame2_b)));
        
        frame2_lab=rgb2lab(frame2);
        frame2_r=squeeze(frame2_lab(:,:,2));
              
        tic;

        tempf1=zeros(size(tempf,1),size(tempf,2));
        tempf1(round(behav.ROI(2)):round(behav.ROI(2)+behav.ROI(4)),round(behav.ROI(1)):round(behav.ROI(1)+behav.ROI(3)))=frame2_r;
        tempf1_s=tempf1;
        thres_red=0.8;
        bwr=tempf1_s>max_framer*thres_red;
        bwr=imdilate(bwr,se);
%         bwr=~bwareaopen(bwr,1000);
        bwrt=bwr;
        
        tempf2=zeros(size(tempf,1),size(tempf,2));
        tempf2(round(behav.ROI(2)):round(behav.ROI(2)+behav.ROI(4)),round(behav.ROI(1)):round(behav.ROI(1)+behav.ROI(3)))=frame2_b;
%         tempf2_s=imgaussfilt(tempf2,3,'FilterSize',[15 15]);
        tempf2_s=tempf2;
        thres_blue=0.8;
        bwb=tempf2_s>max_frameb*thres_blue;
        bwb=imdilate(bwb,se);
%         bwb=~bwareaopen(bwb,1000);
        
        if r_b_diff>10
            bwb(bwrt)=0;
        else
            bwb=bwr;
        end
        bwb=imclose(bwb,se);

        bw=bwr;
        temp1 = rgb2hsv(tempf);
        temp2 = tempf;
        H = temp1(:,:,1);
        S = temp1(:,:,2);
        V = temp1(:,:,3);
        H = mean(mean(H(bw)));
        S = mean(mean(S(bw)));
        V = mean(mean(V(bw)));
        R = temp2(:,:,1);
        G = temp2(:,:,2);
        B = temp2(:,:,3);
        R = mean(mean(R(bw)));
        G = mean(mean(G(bw)));
        B = mean(mean(B(bw)));
        hsvLevelsr = [H+[-.2 .2] S+[-.2 .2] V+[-.2 .2]];
        rgbLevelsr = [R+[-.2 .2] G+[-.2 .2] B+[-.2 .2]];
        statss=regionprops(bw,'centroid');
        centroid_all=reshape([statss.Centroid],2,size([statss.Centroid],2)/2)';
        rcentroid=mean(centroid_all,1);
%         if ~isnan(sum(rcen_p))&&~isempty(sum(rcen_p))&&~isempty(centroid_all)            
%             dis_to_pervious_centroid=sum((centroid_all-repmat(rcen_p,size(centroid_all,1),1)).^2);
%             smallest=find( dis_to_pervious_centroid==min( dis_to_pervious_centroid));
%             centroid_all=[centroid_all;repmat(centroid_all(smallest,:),3,1)];% add weight to the closest one to the point in pervious frame
%             rcentroid=mean(centroid_all,1);
%         else
%             if ~isempty(centroid_all)
%                 rcentroid=mean(centroid_all,1);
%             else
%                 rcentroid=[nan nan];
%             end
%         end
        
       
        bw=bwb;
        temp1 = rgb2hsv(tempf);
        temp2 = tempf;
        H = temp1(:,:,1);
        S = temp1(:,:,2);
        V = temp1(:,:,3);
        H = mean(mean(H(bw)));
        S = mean(mean(S(bw)));
        V = mean(mean(V(bw)));
        R = temp2(:,:,1);
        G = temp2(:,:,2);
        B = temp2(:,:,3);
        R = mean(mean(R(bw)));
        G = mean(mean(G(bw)));
        B = mean(mean(B(bw)));
        hsvLevelsb = [H+[-.2 .2] S+[-.2 .2] V+[-.2 .2]];
        rgbLevelsb = [R+[-.2 .2] G+[-.2 .2] B+[-.2 .2]]; 
        statss=regionprops(bw,'centroid');
        centroid_all=reshape([statss.Centroid],2,size([statss.Centroid],2)/2)';        
        bcentroid=mean(centroid_all,1);
%         if ~isnan(sum(bcen_p))&&~isempty(sum(bcen_p))&&~isempty(centroid_all) 
%             dis_to_pervious_centroid=sum((centroid_all-repmat(bcen_p,size(centroid_all,1),1)).^2);
%             smallest=find(dis_to_pervious_centroid==min(dis_to_pervious_centroid));
%             centroid_all=[centroid_all;repmat(centroid_all(smallest,:),3,1)];% add weight to the closest one to the point in pervious frame
%             bcentroid=mean(centroid_all,1);
%         else
%             if ~isempty(centroid_all)
%                 bcentroid=mean(centroid_all,1);
%             else
%                 bcentroid=[nan nan];
%             end
%         end
        