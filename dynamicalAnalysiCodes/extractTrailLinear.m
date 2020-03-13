function [neuron,idxLRi,idxRLi,idxTrail] = extractTrailLinear(neuron, behav)
    %% Break position into trialNums

    pos = behav.position(:,1);
    neuron.pos = interp1(behav.time, pos, neuron.time);
    
    minP = min(neuron.pos);
    maxP = max(neuron.pos);
    
    edge1 = behav.trackLength/10; %20 - cuts off 10% (was 5%)on each side 
    edge2 = behav.trackLength - behav.trackLength/10; %20
    
    tempPos = neuron.pos;
    tempPos(isnan(tempPos))=-100;
    %%%% the following comments are added by Suoqin on 10-08-2016
    [t1,p1,idx1] = polyxpoly(neuron.time,tempPos,[1 neuron.time(end)],edge1*[1 1]); % find the intersections of the position trace and the boundary. 
    % that is , when the mouse is at left. idx1 is the index of time of the insection, so it records the time for mouse moving to the left
    % idx2 records the time for mouse moving to the right
    idx1 = idx1(:,1);
    [t2,p2,idx2] = polyxpoly(neuron.time,tempPos,[1 neuron.time(end)],edge2*[1 1]);
    idx2 = idx2(:,1);
 
    edgeCrossings = [];
    edgeCrossings(:,1) = [ones(1,length(idx1)) 2*ones(1,length(idx2))];
    edgeCrossings(:,2) = [idx1' idx2'];% in the first column of this matrix, 1 indicates the left, and 2 is the right. 
    edgeCrossings = sortrows(edgeCrossings,2);% sort the matrix based on the second column. so the sorted matrix records the tracjective of movement, 1-2-2-1 means from the left to the right, stay in the right , and then return to the left.
    diffEC = [diff(edgeCrossings(:,1)); 0];
% diffEC=1 means moving from left to the right, -1 means from the right to the
% lift. 
    idx1 = find(diffEC==1);% 
    idx2 = find(diffEC==-1);
    temp = [];
    temp(:,1) = edgeCrossings(sort([idx1' idx2']),2);% the first row indicates the first time for moving from the left to the right. the content is the index time when it moves 
    temp(:,2) = edgeCrossings(sort([idx1' idx2'])+1,2); % the second row indicates the first time fot moving back. the third ros indicates the second time for moving from the left to right
    neuron.trialNum = zeros(neuron.num2read,1);
    for i=1:length(temp)
        neuron.trialNum(temp(i,1):temp(i,2)) = i;% Break position into trialNums, a total of #length(temp) trialNums. i indicates the i-th movement from one side to another
    end
%% extract the time points for each trail
% idxLR = mod(neuron.trialNum,2)==1; % the index for moving from the left to the right
% idxRL = mod(neuron.trialNum+1,2)==1 & neuron.trialNum~=0;% for moving from the right to the left
idxLRi = [];idxRLi = [];
for trialNum = 1:2:max(neuron.trialNum)
    iidx = (neuron.trialNum == trialNum);
    idxLRi = [idxLRi, iidx];
end
for trialNum = 2:2:max(neuron.trialNum)
    iidx = (neuron.trialNum == trialNum);
    idxRLi = [idxRLi, iidx];
end
idxTrail = [idxLRi,idxRLi];
idxTrail = logical(idxTrail);  
end


