function displayIndividualTrails(neuron0,segDisplay,thresh)
%% display the trail and overlay the neuron activity onto trajectories
figure
% num = 6;  % only dislay the first six trails in each direction
k = 0;
for segNum = segDisplay%size(neuron1.S,1)
    k = k+1;
    idx = neuron0.S(segNum,:) > thresh(segNum); idx = idx';
    subplot(length(segDisplay),2,2*k-1)
    for trialNum = 1:2:max(neuron0.trialNum)
        iidx = (neuron0.trialNum == trialNum);
        plot(neuron0.time(iidx),neuron0.pos(iidx ),'k','linewidth',1.5);
        hold on
        plot(neuron0.time(idx&iidx),neuron0.pos(idx&iidx),'r.')
    end
    if segNum ~= max(segDisplay)
        set(gca,'Xtick',[])
    end
    title(['Neuron ' num2str(segNum) '. Direction 1'],'FontSize',10)
    hold off
    
    subplot(length(segDisplay),2,2*k)
    for trialNum = 2:2:max(neuron0.trialNum)
        iidx = (neuron0.trialNum == trialNum);
        plot(neuron0.time(iidx),neuron0.pos(iidx),'k','linewidth',1.5);
        hold on
        plot(neuron0.time(idx&iidx),neuron0.pos(idx&iidx),'r.')
    end
    if segNum ~= max(segDisplay)
        set(gca,'Xtick',[])
    end
    title(['Neuron ' num2str(segNum) '. Direction 2'],'FontSize',10)
    hold off
end