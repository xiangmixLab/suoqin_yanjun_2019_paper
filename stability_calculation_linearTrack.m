neuronIndividuals = neuronIndividualsf;behavIndividuals = behavIndividualsf;
%% ensemIfo- ensemble information// ensemIfo1 condition1; RL1 right to left condition 1
binSize = 10;
ensemIfo = cell(1,length(neuronIndividuals));
ensemIfoLR = ensemIfo;ensemIfoRL = ensemIfo;ensemIfoLRodd = ensemIfo;ensemIfoRLodd = ensemIfo;
ensemIfoLReven = ensemIfo;ensemIfoRLeven = ensemIfo;
ensemIfoHalfFirst = cell(1,length(neuronIndividuals)); ensemIfoHalfSecond = cell(1,length(neuronIndividuals));
for i = 1:length(neuronIndividuals)
    [ensemIfo{i},ensemIfoLR{i},ensemIfoRL{i},ensemIfoLRodd{i},ensemIfoRLodd{i},ensemIfoLReven{i},ensemIfoRLeven{i}] = ...
        calculatingEnsembleActivityLinearTrack2(neuronIndividuals{i},behavIndividuals{i},thresh,'S',binSize);
    [ensemIfoHalfFirst{i},ensemIfoHalfSecond{i}] = calculatingEnsembleActivityLinearTrackHalfFirstSecond(neuronIndividuals{i},behavIndividuals{i},thresh,'S',binSize);
end

%% calculate the correlation between odd and even trails
% corr_diag_ave = zeros(1,length(ensemIfoLRodd));
corr_diag_ave = zeros(size(ensemIfoRLodd{1}.FR,2),length(ensemIfoLRodd));
corrIndividuals = cell(1,length(ensemIfoLRodd));
for i = 1:length(ensemIfoLRodd)
    RLodd = ensemIfoLRodd{i};RLeven = ensemIfoLReven{i}; % trails from left to right
    % RLodd.FR(isnan(RLodd.FR))=0;RLeven.FR(isnan(RLeven.FR))=0;
    d = pdist2(RLodd.FR',RLeven.FR','correlation');
    corr = 1-d;corr(isnan(corr))=0;
%     corr_diag_ave(i) = mean(diag(corr));
    corr_diag_ave(:,i) = diag(corr);
    corrIndividuals{i} = corr;
end
corr_diag_LR = corr_diag_ave;corrIndividuals_LR = corrIndividuals;

corr_diag_ave = zeros(size(ensemIfoRLodd{1}.FR,2),length(ensemIfoLRodd));
corrIndividuals = cell(1,length(ensemIfoLRodd));
for i = 1:length(ensemIfoLRodd)
    RLodd = ensemIfoRLodd{i};RLeven = ensemIfoRLeven{i}; % trails from left to right
    % RLodd.FR(isnan(RLodd.FR))=0;RLeven.FR(isnan(RLeven.FR))=0;
    d = pdist2(RLodd.FR',RLeven.FR','correlation');
    corr = 1-d;corr(isnan(corr))=0;
%     corr_diag_ave(i) = mean(diag(corr));
    corr_diag_ave(:,i) = diag(corr);
    corrIndividuals{i} = corr;
end
corr_diag_RL = corr_diag_ave;corrIndividuals_RL = corrIndividuals;


%save corr_odd_even_M2016.mat corrIndividuals_LR corrIndividuals_RL corr_diag_LR corr_diag_RL 


corr_half_diag_ave = zeros(size(ensemIfoHalfFirst{1}.FR,2),length(ensemIfoHalfFirst));
corr_halfIndividuals = cell(1,length(ensemIfoHalfFirst));
for k = 1:length(neuronIndividuals)
%     FR_first = ensemIfoHalfFirst{k}.FR(idx_PC,:); % firing rate from left to right
%     FR_second = ensemIfoHalfSecond{k}.FR(idx_PC,:); %firing rate from left to right
    FR_first = ensemIfoHalfFirst{k}.FR; % firing rate from left to right
    FR_second = ensemIfoHalfSecond{k}.FR; %firing rate from left to right
    d = pdist2(FR_first',FR_second','correlation');
    corr0 = 1-d;corr0(isnan(corr0))=0;
%     corr_diag_ave(i) = mean(diag(corr));
    corr_half_diag_ave(:,k) = diag(corr0);
    corr_halfIndividuals{k} = corr0;
end
%save spatial_corr_half_first_second_M2016_placeCells_R1.mat corr_half_diag_ave corr_halfIndividuals sessionsSlected

save stability_analysis_R2.mat corrIndividuals_LR corrIndividuals_RL corr_diag_LR corr_diag_RL corr_half_diag_ave corr_halfIndividuals


