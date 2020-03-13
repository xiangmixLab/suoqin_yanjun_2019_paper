function thresh = determiningFiringEventThresh(neuron,temp)
% maxS = max(neuron.S,[],2);
% figure;hist(maxS)
% % thresh = mean(maxS)+2*std(maxS);
% thresh = mean(maxS)-1*std(maxS);
% hold on
% line([thresh thresh],[0 30],'Color','red','LineStyle','--')
% xlabel('Peak of spikes','FontSize',10)
% ylabel('Count','FontSize',10)

% figure;
% for i = 1:length(neuronIndividuals)
%     subplot(1,length(neuronIndividuals),i)
%     maxS = max(neuronIndividuals{i}.S,[],2);
%     hist(maxS)
% 
% hold on
% line([maxS maxS],[0 30],'Color','red','LineStyle','--')
% xlabel('Peak of spikes','FontSize',10)
% ylabel('Count','FontSize',10)
% end
if ~exist('temp','var'), temp = 'S'; end
if strcmpi(temp,'trace')
    data = neuron.trace;
elseif strcmpi(temp,'S')
    data = neuron.S;
end
maxS = max(data,[],2);
% thresh = 0.1*maxS; % a vector
thresh = 0.1*maxS; % a vector