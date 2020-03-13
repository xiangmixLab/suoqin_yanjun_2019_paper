function [D12, D32, P12, P32] = distributionTest_bootstrap(infoScorebootAll,sessionsSlected)

infoScore = cell(1,length(sessionsSlected));
for sessionIndex = 1:length(sessionsSlected)
    infoScore{sessionIndex} = infoScorebootAll{1,sessionsSlected(sessionIndex)}(:,2:end);
end

D12 = infoScore{1} - infoScore{2};
D32 = infoScore{3} - infoScore{2};

Q12 = quantile(D12,[0.025 0.975],2);P12 = Q12(:,1).*Q12(:,2) > 0; 
Q32 = quantile(D32,[0.025 0.975],2);P32 = Q32(:,1).*Q32(:,2) > 0; 
