function generateAgeSubsampleAndCheckPrecision( )

    
    fileName = 'randomSizeAgePrecision.mat';
    
    load('easyFormatHumanData_withAges.mat');
    load('subset.mat');
    load('humanOntologyObject.mat');
    
    sizesOfClasses = sum(ages,1);
    numberOfPeople = size(experimentsSubjectMatrixLogical,2);
    numberOfRepetitions = 10000;
    numberOfEra = size(agesDescription,1);

    experimentsDataMatrixNormalized = normalizeMatrixByDim(experimentsDataMatrix, 1);
    experimentsDataMatrix =experimentsDataMatrix;
    experimentsSubjectMatrixLogical = experimentsSubjectMatrixLogical;
    
    regionSampleTreshold = 50;
    
    regionColors = humanOntology.structureColors(subsetNodesIndex,:);
    regionNames = humanOntology.structureLabels(subsetNodesIndex,4);
    experimentRegion = getRegionOfGrossStructures(subsetNodesIndex, experimentsLocationMatrix, humanOntology);

    experimentsInRegion = sum(experimentRegion,1);
    aboveTreshold = (experimentsInRegion > regionSampleTreshold);
    experimentRegion = experimentRegion(:,aboveTreshold);
    regionNames = regionNames(aboveTreshold);
    regionColors = regionColors(aboveTreshold,:);
    experimentsInRegion = experimentsInRegion(aboveTreshold);

    numberOfRegions = size(experimentRegion,2);
    
    
    
    [categoriesScores, distanceBetweenGroups] = calcRegionBetweenSubjectDist(experimentsSubjectMatrixLogical, experimentsDataMatrixNormalized, experimentRegion);
    
    percentInCategory = nan(numberOfRepetitions,numberOfPeople, numberOfRegions,length(sizesOfClasses) );
    percentInCategoryUsingGenes = nan(numberOfRepetitions,numberOfPeople, numberOfRegions,length(sizesOfClasses) );
    
    parfor i = 1:length(sizesOfClasses)
        sizeOfClass = sizesOfClasses(i);
        [percentInCategory(:,:,:,i) , percentInCategoryUsingGenes(:,:,:,i)] = getClusteringScoresOfArandomSet(numberOfRepetitions, sizeOfClass, distanceBetweenGroups, experimentsDataMatrix, experimentsSubjectMatrixLogical);
%         printPercentCounter(i, length(sizesOfClasses));
    end
    
    save(fileName,'percentInCategory','percentInCategoryUsingGenes','sizesOfClasses','experimentsDataMatrix', 'experimentsSubjectMatrixLogical','experimentRegion','ages','agesDescription','selectedProbesData','regionNames','regionColors');
    
end

function [percentInCategory,percentInCategoryWithGenes] = getClusteringScoresOfArandomSet(numberOfRepetitions,sizeOfClass, distanceBetweenGroups, experimentsDataMatrix, experimentsSubjectMatrixLogical)
    totalNumberOfGenes = size(experimentsDataMatrix,2);
    numberOfRegions = size(distanceBetweenGroups,1);
    numberOfCategories = size(experimentsSubjectMatrixLogical,2);
    
    percentInCategory = nan(numberOfRepetitions,numberOfCategories,numberOfRegions);
    percentInCategoryWithGenes = nan(numberOfRepetitions,numberOfCategories,numberOfRegions);
    %l2PercentInCategory = nan(numberOfRepetitions,1);
%    fprintf('\nj completed    ');
    randomPermMatrix = false(numberOfRepetitions, totalNumberOfGenes);
    for  i=1:numberOfRepetitions
       randomSubSet = randperm(totalNumberOfGenes,sizeOfClass);
       randomPermMatrix (i,randomSubSet) = true;
    end
    
    parfor m =1:numberOfRegions
        currentAreaDistanceBetweenGroups = squeeze(distanceBetweenGroups(m,:,:,:));
        distanctUsingAllFeatures = sum(currentAreaDistanceBetweenGroups,3);
        distanctUsingAllFeatures = squeeze( distanctUsingAllFeatures );
    
        for j=1:numberOfRepetitions
            randomSubSet = randomPermMatrix(j,:);
            % first sum over all genes
            
            %distancesWithoutTheSubset = currentAreaDistanceBetweenGroups(:,:,~randomSubSet);
            %percentInCategory(j,:,m) = clusteringScore( distancesWithoutTheSubset );
            
            distancesWithOnlyTheSubset = currentAreaDistanceBetweenGroups(:,:,randomSubSet);
            [percentInCategory(j,:,m) ,percentInCategoryWithGenes(j,:,m)] = clusteringScoreAlternatrive(distancesWithOnlyTheSubset,distanctUsingAllFeatures);
        end
    end
end

function percentInCategory = getScoresOfArandomSet(numberOfRepetitions,sizeOfClass, expressionDistances, experimentsDataMatrix, experimentsSubjectMatrixLogical)
    totalNumberOfGenes = size(experimentsDataMatrix,2);
    numberOfCategories = size(experimentsSubjectMatrixLogical,2);
    
    percentInCategory = nan(numberOfRepetitions,numberOfCategories);
    %l2PercentInCategory = nan(numberOfRepetitions,1);
%    fprintf('\nj completed    ');
    parfor j=1:numberOfRepetitions
        randomSubSet = randperm(totalNumberOfGenes,sizeOfClass);
        newDistanceMatrix = removeFeaturesFromDistance(randomSubSet, expressionDistances, experimentsDataMatrix);

        %newDistanceMatrix = squareform(newDistanceMatrix,'tomatrix');

        percentInCategory(j,:) = calcPrecisionAt1( newDistanceMatrix, experimentsSubjectMatrixLogical, true);
 %       printPercentCounter(j, numberOfRepetitions );
    end
end

function newDistanceMatrix = removeFeaturesFromDistance(featureToRemoveIndices, distanceMatrix, samplesInRn)
    numberOfSamples = size(samplesInRn,1);
    numberOfFeatures = size(samplesInRn,2);
    
    onlyTheRemovedFeatures = samplesInRn(:, featureToRemoveIndices);
    %removeFeaturesDistance = pdist2(onlyTheRemovedFeatures,onlyTheRemovedFeatures,'euclidean');
    %removeFeaturesDistance = distancesUsingPdist(onlyTheRemovedFeatures);

    squaredDistances = unSquaredPdist(onlyTheRemovedFeatures);
    
    newDistanceMatrix = distanceMatrix.^2 - squaredDistances;

end

function distances = distancesUsingPdist(x)
        distanceVector = pdist(x,'euclidean');
        distances = squareform(distanceVector);
end
function squaredDistances = unSquaredPdist(x)
    Qx=repmat(dot(x,x,2),1,size(x,1));
    squaredDistances = (Qx+Qx'-2*(x*x'));
end