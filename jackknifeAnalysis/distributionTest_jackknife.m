function [Q, P] = distributionTest_jackknife(infoScorebootAll,sessionsSelected)
% sessionsSelected is the index of three sessions: Ctrl, CNO, and Pctrl
mu = [];Var = [];
for sessionIndex = 1:length(sessionsSelected)
    infoScore = infoScorebootAll{1,sessionsSelected(sessionIndex)}(:,2:end); % get all the info score from randomly sampling
    mu(:,sessionIndex) = mean(infoScore,2);
    Var(:,sessionIndex) = var(infoScore,1,2)*(size(infoScore,2)-1);
end
% construct the jackknife-based t-test statistics
Q12 = (mu(:,1) - mu(:,2))./sqrt(Var(:,1) + Var(:,2));
Q32 = (mu(:,3) - mu(:,2))./sqrt(Var(:,3) + Var(:,2));
Q = [Q12, Q32];
P = normcdf(abs(Q),'upper')*2;