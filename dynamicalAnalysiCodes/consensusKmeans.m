function CM = consensusKmeans(neuron0,threshCluster,K,N)
% function CM = consensusKmeans(idxLRi,idxRLi,neuron0,threshCluster,K,N)
%% performing clustering using kmeans+consensus clustering method
% using neuron.trace
% K = 5; % initial guess of the number of clusters
% N = 1000; 
% n = size(idxLRi,2)+size(idxRLi,2);
rng default  % For reproducibility
% idxTrail = [idxLRi,idxRLi];
% idxTrail = logical(idxTrail);
CM = zeros(size(neuron0.trace,1));
for i = 1:N
%     trialNum = randperm(n,1);
%     datai = neuron0.trace(:,idxTrail(:,trialNum));
    datai = neuron0.trace(:,logical(randi([0 1],size(neuron0.trace,2),1)));
    idxD = find(max(datai,[],2) < threshCluster);
    idxC = setdiff(1:size(datai,1),idxD);
    datai(idxD,:) = [];groupD = [idxD,(K+1)*ones(length(idxD),1)];
    try
    group = kmeans(datai,K,'dist','corr','rep',10,'disp','final');
    groupC = [idxC(:),group];
    group = [groupC; groupD];
    [~,idx] = sort(group(:,1)); group = group(idx,2);
    for k = 1:length(unique(group))
        comb = nchoosek(find(group == k),2);
        if length(comb) == 1
            continue;
        end
        linearInd = sub2ind(size(CM), comb(:,1),comb(:,2));
        CM(linearInd) = CM(linearInd) + 1;
    end
    catch
        continue;
    end
end
CM = max(CM,CM')/N;

cgo = clustergram(CM,'Standardize',3,'Linkage','complete');
set(cgo,'Colormap',redbluecmap);