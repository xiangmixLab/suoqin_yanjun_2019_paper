function [data1_mu,data1_sd, data2_mu,data2_sd, data3_mu,data3_sd, data4_mu,data4_sd] = firingRateCalculation_object(sumFiringRateObjectAllAll_saline,sumFiringRateObjectAllAll_CNO)

data1 = [];
for i = 1:size(sumFiringRateObjectAllAll_saline,1)
    data1 = [data1;sumFiringRateObjectAllAll_saline{i,2}];
end
data1_mu = mean(data1); data1_sd = std(data1); 
data2 = [];
for i = 1:size(sumFiringRateObjectAllAll_CNO,1)
    data2 = [data2;sumFiringRateObjectAllAll_CNO{i,2}];
end
data2_mu = mean(data2); data2_sd = std(data2); 
data3 = [];
for i = 1:size(sumFiringRateObjectAllAll_saline,1)
    data3 = [data3;sumFiringRateObjectAllAll_saline{i,3}];
end
data3_mu = mean(data3); data3_sd = std(data3); 
data4 = [];
for i = 1:size(sumFiringRateObjectAllAll_CNO,1)
    data4 = [data4;sumFiringRateObjectAllAll_CNO{i,3}];
end
data4_mu = mean(data4); data4_sd = std(data4); 

