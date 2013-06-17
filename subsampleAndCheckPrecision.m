function subsampleAndCheckPrecision(sizesOfClasses, fileName)

    if ~exist('fileName','var')
        fileName = 'randomSizePrecision.mat';
    end

    geneGo = load('humanGene2GoMatrix.mat', 'go_genes_mat', 'cat_ids', 'geneNames', 'aspects');
    load('expressionRawDataWithOntology.mat', 'selectedProbesData');

    if exist('sizesOfClasses','var')
        sizesOfClasses = setdiff(sizesOfClasses, unique(sum(geneGo.go_genes_mat,2)) );
    else
        sizesOfClasses = unique(sum(geneGo.go_genes_mat,2));
    end
    
    load('dataWithDistance.mat');
    
    
    numberOfPeople = 6;
    numberOfRepetitions = 1000;

    experimentsSubjectMatrixLogical = false(size(experimentsSubjectMatrix,1), numberOfPeople);
    for i =1:numberOfPeople
        experimentsSubjectMatrixLogical(:,i) = experimentsSubjectMatrix == i;
    end
    
    
    
    
    %find genes which appear both in GO and in Allen
    [commonGenes,commonGeneInGo,commonGeneInAllen] = intersect(geneGo.geneNames, selectedProbesData.gene_symbols);
    go_gene_mat = geneGo.go_genes_mat(:, commonGeneInGo);
    geneNames = geneGo.geneNames(commonGeneInGo);
    go_cat_without_genes = sum(go_gene_mat,2) == 0;
    cat_ids = geneGo.cat_ids(~go_cat_without_genes);
    go_gene_mat = go_gene_mat(~go_cat_without_genes, :);
    experimentsDataMatrix = experimentsDataMatrix(:, commonGeneInAllen);
    aspects = geneGo.aspects(~go_cat_without_genes);
    
    expressionDistances = sqrt(    unSquaredPdist(experimentsDataMatrix)    );
    sizesOfClasses = unique(sum(go_gene_mat,2));
    numberOfSamples = size(experimentsSubjectMatrixLogical,1);
    numberOfCategories = size(experimentsSubjectMatrixLogical,2);
    totalNumberOfGenes = size(experimentsDataMatrix,2);
    numberOfRegions = size(experimentRegion,2);
    assert(  numberOfSamples == size(experimentsDataMatrix,1) );

    % first remove randomly sampled genes
%     meanScores = nan(length(sizesOfClasses),1);
%     meanL2Scores = nan(length(sizesOfClasses),1);
%     stdScores = nan(length(sizesOfClasses),1);
%     stdL2Scores = nan(length(sizesOfClasses),1);
%     
    
    diagIndices = 1:size(expressionDistances,1)+1:numel(expressionDistances); 
    expressionDistances(diagIndices) = inf;
    
    %distanceAsAVector = squareform(expressionDistances,'tovector');

    fprintf('\ni completed    ');
    
    
    [categoriesScores, distanceBetweenGroups] = calcRegionBetweenSubjectDist(experimentsSubjectMatrixLogical, experimentsDataMatrix, experimentRegion);
    
    percentInCategory = nan(numberOfRepetitions,numberOfCategories, numberOfRegions,length(sizesOfClasses) );
    percentInCategoryUsingGenes = nan(numberOfRepetitions,numberOfCategories, numberOfRegions,length(sizesOfClasses) );
    
    size(distanceBetweenGroups)
    for i = 1:length(sizesOfClasses)
        tic
        sizeOfClass = sizesOfClasses(i);

        [percentInCategory(:,:,:,i) , percentInCategoryUsingGenes(:,:,:,i)] = getClusteringScoresOfArandomSet(numberOfRepetitions, sizeOfClass, distanceBetweenGroups, experimentsDataMatrix, experimentsSubjectMatrixLogical);
%         toc;
       % fprintf('     ');
        printPercentCounter(i, length(sizesOfClasses));
       % fprintf('\n');
    end
    
    save(fileName,'percentInCategory','percentInCategoryUsingGenes','sizesOfClasses','experimentsDataMatrix', 'experimentsSubjectMatrixLogical','experimentRegion','cat_ids', 'aspects','go_gene_mat','geneNames');
    
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
    
    for m =1:numberOfRegions
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