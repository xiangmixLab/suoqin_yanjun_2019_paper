function discriminationIndexComparison_differentRadius(r, DIall_mu_saline,DIall_sd_saline, DIall_mu_CNO, DIall_sd_CNO,sessionID,sessionName)
colors_conditions = [228,26,28;55,126,184]/255;
hFig = figure('position', [200, 200, 190,170]);
h(1) = errorbar(r,DIall_mu_saline(:,sessionID),DIall_sd_saline(:,sessionID));
hold on
h(2) = errorbar(r,DIall_mu_CNO(:,sessionID),DIall_sd_CNO(:,sessionID));
for i = 1:2
    h(i).Color = colors_conditions(i,:);
    h(i).LineWidth = 1;
end
box off
legend({'Saline', 'CNO'},'FontName','Arial')
xlim([1.8 5.2])
set(gca,'xtick',r)
set(gca,'Xticklabel',{'1.5','3','4.5','6'})
xlabel('Radius (cm)','FontName','Arial')
ylabel({'Neural discrimination', 'index (event rate,%)'},'FontName','Arial')
set(gca,'FontSize',12,'FontName','Arial')
title(sessionName,'FontSize',12,'FontName','Arial')