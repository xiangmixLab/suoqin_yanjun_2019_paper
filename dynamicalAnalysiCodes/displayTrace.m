function displayTrace(neuron0,trialNum, group,idxTrail,colorClusters)
%% display the trace of each cluster
% trialNum = [1:10]; % the trails for displaying, e.g. the first ten trails
figure
optimalK = length(unique(group));
ha = tight_subplot(optimalK,1,[.03 .03],[.05 .05],[.08 .01]);
for i = 1:optimalK
    datai = [];
    for j = 1:length(trialNum)
        datai = [datai,neuron0.trace(:,idxTrail(:,trialNum(j)))];
    end
    axes(ha(i))
    %     subplot(optimalK,2,2*i-1);
    plot(1:size(datai,2),datai((group == i),:));
%     plot(1:size(datai,2),mean(datai((group == i),:)),'k');
    title(['Num of neurons:',num2str(length(find(group==i)))],'FontSize',8)
    set(gca,'FontSize',8)
    axis tight
    if i ~= optimalK, set(gca,'Xtick',[]);end
end

figure
optimalK = length(unique(group));
ha = tight_subplot(optimalK,1,[.03 .03],[.05 .05],[.08 .01]);
for i = 1:optimalK
    datai = [];
    for j = 1:length(trialNum)
        datai = [datai,neuron0.trace(:,idxTrail(:,trialNum(j)))];
    end
    axes(ha(i))
    %     subplot(optimalK,2,2*i-1);
    idx = find(group == i);
    [~,order] = sort(max(datai(idx,:),[],2),'ascend');
    idx = idx(order);
    for j = 1:nnz(group == i)
        increment = j*10;
        plot(1:size(datai,2),datai(idx(j),:)+increment,'color',colorClusters(i,:),'LineWidth',0.5);
        hold on
    end
%     plot(1:size(datai,2),mean(datai((group == i),:)),'k');
%     title(['Num of neurons:',num2str(length(find(group==i)))],'FontSize',8)
    set(gca,'FontSize',8)
    axis tight
%     if i ~= optimalK, set(gca,'Xtick',[]);end
    axis off
end

figure
optimalK = length(unique(group));
% groupCenter = zeros(optimalK,size(smoothCurve,2));
ha = tight_subplot(optimalK,1,[.03 .03],[.05 .05],[.08 .01]);
for i = 1:optimalK
    datai = [];
    for j = 1:length(trialNum)
        datai = [datai,neuron0.trace(:,idxTrail(:,trialNum(j)))];
    end
    axes(ha(i))
    %     subplot(optimalK,2,2*i-1);
%     plot(1:size(datai,2),datai((group == i),:),'k');
        Q = quantile(datai((group == i),:),[0.25, 0.5, 0.75]);
        groupCenter(i,:) = 0.25*Q(1,:)+0.5*Q(2,:)+0.25*Q(3,:);
%     plot(1:size(datai,2),groupCenter(i,:),'color',colorClusters(i,:),'LineWidth',1);
    plot(1:size(datai,2),mean(datai((group == i),:)),'color',colorClusters(i,:),'LineWidth',1);
%     title(['Num of neurons:',num2str(length(find(group==i)))],'FontSize',8)
    set(gca,'FontSize',8)
    axis tight
    ylim([0 3.5])
    axis off
    if i ~= optimalK, set(gca,'Xtick',[]);end
end
