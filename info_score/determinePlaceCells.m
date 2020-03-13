function [place_cells,infoScoreThresh] = determinePlaceCells(infoScorenull,infoScoreboot,neuron_lowFR)
if ~exist('neuron_lowFR','var') || isempty(neuron_lowFR)
    neuron_lowFR = [];
end
infoScore = infoScoreboot;

infoScore(neuron_lowFR,:) = [];

infoScoreThresh = quantile(infoScore,0.95,2);

TinfoScorenull = table([1:length(infoScorenull)]',infoScorenull,infoScoreThresh,'VariableNames',{'neuron','infoPerSecond','thresh'});
TinfoScorenull(neuron_lowFR,:) = [];
TinfoScorenull2 = sortrows(TinfoScorenull,{'infoPerSecond'},{'descend'});
place_cells = TinfoScorenull2.neuron(TinfoScorenull2.infoPerSecond > TinfoScorenull2.thresh);
