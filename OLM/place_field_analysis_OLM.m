function place_field_analysis_OLM(filefolder,mice)

resultsfile = fullfile(filefolder,[mice,'_results_bin15_3SD.mat'])
load(resultsfile, 'FiringRateAll', 'posObjectsAll')

% load place cells identified using individual shuffling
load(fullfile(filefolder,'results',['place_cells_',mice,'_individualShuffle.mat']),'place_cellsAll')

nSessions = size(FiringRateAll,2);
radius = 3 - 1; % event rate associated with the object was quantified as the average event rate across bins within (radius-1)bins surrounding the object.

PFresponseAll = cell(1, nSessions); peakERresponseAll = cell(1, nSessions);
peakRateAll = cell(1, nSessions);
for sessionIndex = 1:nSessions
    firingRate = FiringRateAll(:,sessionIndex);
    %     segments = place_cellsAll{sessionIndex}; % use place cells defined in each session
    segments = place_cellsAll{length(place_cellsAll)-1}; % use place cells defined in training
    
    posObjects = posObjectsAll{sessionIndex};
    idx_posObjects = zeros(size(posObjects,1), 2);
    if sum(posObjects(:))~=0
        for i = 1:size(posObjects,1)
            idx_posObjects(i,:) = [size(firingRate{1},1)-posObjects(i,2)+1, posObjects(i,1)];
        end
    else
        continue;
    end
    
    PFresponse = zeros(length(segments), 3);
    peakERresponse = zeros(length(segments), 3);
    % the first and second column indicates the response to object 1 and 2,respectively.
    % the third column indicates no response to the two objects,but with firing field
    peakRate = zeros(length(segments),1);
    
    for i = 1:length(segments)
        rate_mat0 = firingRate{segments(i)}; % rate map of individual neuron segments(i)
        if isempty(rate_mat0), continue; end
        peakRate(i) = max(rate_mat0(:));
        rate_mat = filter2DMatrices(rate_mat0, 1); % smooth the rate map
        % Create AutoCorr
        autocorr = Cross_Correlation(rate_mat, rate_mat);
        
        %Find AutoMaxInds
        auto_max_inds = FindAutoMaxInds(autocorr);
        if isnan(auto_max_inds), continue; end
        
        % Find Place field radius
        % Finds size of fields (PF_radius) using autocorrelation by taking around 70% of the half distance from the autocorr center to the closest local maximum.
        PF_radius = findPlaceFieldRadius(autocorr, auto_max_inds);
        
        % Finds local maxima (max_inds) of smoothed rate map
        max_inds= FindMaxIndsRateMap(rate_mat);
        strength= 1.9;
        max_inds= RemoveTooCloseMaxInds(max_inds, PF_radius, rate_mat, strength); % removes any closer than the PF_radius to remove points too close together
        
        max_inds_delete = []; % delete the fields with peak rate less than 25% of maximum of rate map
        for ii = 1:size(max_inds)
            if rate_mat(max_inds(ii,1), max_inds(ii,2)) < 0.25*max(rate_mat(:))
                max_inds_delete = [max_inds_delete,ii];
            end
        end
        max_inds(max_inds_delete,:) = []; %  the coordinates of peak of each field (each row)
        
        D = pdist2(idx_posObjects,max_inds);[idx, idy] = find(D <= radius*sqrt(2));
        if ~isempty(idx)
            for ii = 1:length(idx)
                PFresponse(i,idx(ii)) = 1;
            end
            max_inds = max_inds(idy,:);
            % Find peak firing rate at max_inds
            peak_rates = findPeakRates(max_inds,rate_mat);
            
            if length(unique(idx)) == 1
                peakERresponse(i, unique(idx)) = grpstats(peak_rates, idx);
            else
                peakERresponse(i, 1:2) = grpstats(peak_rates, idx);
            end
            
        else
            PFresponse(i,3) = 1;
            peak_rates = findPeakRates(max_inds,rate_mat);
            peakERresponse(i, 3) = mean(peak_rates);
        end
    end
    PFresponseAll{sessionIndex} = PFresponse;
    peakERresponseAll{sessionIndex} = peakERresponse;
    peakRateAll{sessionIndex} = peakRate;
end

% save(fullfile(filefolder,'results',['placeFieldResponse_',mice,'_radius3_placeCells.mat']), 'PFresponseAll','peakERresponseAll')
save(fullfile(filefolder,'results',['placeFieldResponse_',mice,'_radius3_placeCells_basedTraining.mat']), 'PFresponseAll','peakERresponseAll','peakRateAll')
