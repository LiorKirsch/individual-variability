function percentInCategory = calcPrecisionAt1( distanceMatrixBetweenSamples, experimentCategories, dontChangeToInfFlag)


    % the distance between a point and itself is zero so we need to remove this
    
    if ~exist('dontChangeToInfFlag','var')
        diagIndices = 1:size(distanceMatrixBetweenSamples,1)+1:numel(distanceMatrixBetweenSamples); 
        distanceMatrixBetweenSamples(diagIndices) = inf;
    end
    
    
    [~, minIndices] = min(distanceMatrixBetweenSamples,[],1);

    nearestNeighborGroup = experimentCategories(minIndices,:);
    numberOfCategories = size(experimentCategories,2);
    
    samplesWhoseNearestNeighborIsInTheSameCateogry = false(size(experimentCategories));
    for i =1:numberOfCategories
        samplesInRegion = experimentCategories(:,i);
        nearestNeighborInRegion = nearestNeighborGroup(:,i);
        samplesWhoseNearestNeighborIsInTheSameCateogry(:,i) = samplesInRegion & nearestNeighborInRegion;
    end

    percentInCategory = sum(samplesWhoseNearestNeighborIsInTheSameCateogry,1) ./ sum(experimentCategories,1);
end