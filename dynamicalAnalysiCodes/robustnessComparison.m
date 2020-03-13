clc;clear
% N = 50; sigma = [0.95, 0.90, 0.80, 0.75, 0.70];
% N = 50; sigma = [0.80, 0.75, 0.70];
N = 50; sigma = [0.90,0.80, 0.70];
methods = {'scEnergy','monocleICA','monocleDDRTree','tscan','dptnoRoot'};
measure = {'PRS','Pearson'};
% data = zeros(length(sigma),length(methods),N);
% for i = 1:length(methods)
%     load(['robustnessCellSampling_',measure{1},'_',methods{i},'.mat'])
%     if strcmp(methods{i},'tscan'),Kcorr = KcorrPath{1};end
%     data(:,i,:) = abs(Kcorr);
%     clear Kcorr
% end
data = zeros(length(sigma),length(methods),nchoosek(N,2));
for i = 1:length(methods)
    load(['robustnessCellSamplingPairwise_',measure{2},'_',methods{i},'.mat'])
    if strcmp(methods{i},'tscan'),Kcorr = KcorrPath{1};end
    %     Kcorr(:,:,1:2) = [];
    Kcorr(:,:,[1 4]) = [];
    r = zeros(nchoosek(N,2),length(sigma));
    for k = 1:length(sigma)
        r(:,k) = nonzeros(triu(Kcorr(:,:,k),1));
    end
    data(:,i,:) = abs(r');
    clear Kcorr
end
figure
h = boxplot2(data);
% Alter linestyle and color
cmap = get(0, 'defaultaxescolororder');
cmap(2,:) = cmap(1,:);
cmap(1,:) = [1 0 0];
for ii = 1:length(methods)
    structfun(@(x) set(x(ii,:), 'color', cmap(ii,:), ...
        'markeredgecolor', cmap(ii,:)), h);
end
set([h.lwhis h.uwhis], 'linestyle', '-');
set(h.out, 'marker', '.');
set(gca,'FontSize',8)
xlim([0.5 length(sigma)+0.5])
xlabels = cellstr(num2str(sigma'*100,'%.0f%%'));
set(gca,'Xtick',1:length(sigma))
set(gca,'XtickLabel',xlabels,'FontSize',10,'FontName','Arial')
ylabel('Robustness','FontSize',10,'FontName','Arial')
xlabel('Subsampling','FontSize',10,'FontName','Arial')
title('Pseudotime reconstruction score','FontSize',10,'FontName','Arial');
title('Pearson correlation','FontSize',10,'FontName','Arial');

