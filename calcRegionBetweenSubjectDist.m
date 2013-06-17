function [categoriesScores, distanceBetweenGroups] = calcRegionBetweenSubjectDist(experimentsSubjectMatrixLogical, experimentsDataMatrix, experimentRegion)
    numberOfPeople = size(experimentsSubjectMatrixLogical,2);
    numberOfSamples = size(experimentsSubjectMatrixLogical,1);
    numberOfRegions = size(experimentRegion,2);
    numberOfGenes = size(experimentsDataMatrix,2);
    assert(  numberOfSamples == size(experimentsDataMatrix,1) );

    
    categoriesScores = nan(numberOfRegions,numberOfPeople);
    distanceBetweenGroups = nan(numberOfRegions,numberOfPeople,numberOfPeople, numberOfGenes);
    for m =1:numberOfRegions
        currentAreaSamplesIndices = experimentRegion(:,m);
        currentAreaDataMatrix = experimentsDataMatrix(currentAreaSamplesIndices,:);
        currentSubjectMatrix = experimentsSubjectMatrixLogical(currentAreaSamplesIndices,:);
        currentDist = distBetweenAndWithinGroups(currentAreaDataMatrix, currentSubjectMatrix);
        distanceBetweenGroups(m,:,:,:) =  currentDist;
        categoriesScores(m,:) = clusteringScore( currentDist );
    end
end
