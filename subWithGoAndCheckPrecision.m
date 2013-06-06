function subWithGoAndCheckPrecision()



    geneGo = load('humanGene2GoMatrix.mat', 'go_genes_mat', 'cat_ids', 'geneNames', 'aspects');
    
    load('dataWithDistance.mat','experimentsSubjectMatrix', 'experimentsDataMatrix', 'expressionDistances','experimentRegion','regionNames');
    load('expressionRawDataWithOntology.mat', 'selectedProbesData');
    
%     randomSizePrecision =  load('randomSizePrecision.mat','percentInCategory','sizesOfClasses');
%     
%     [meanPerPerson, stdPerPerson] = calcNormalForEachPerson(randomSizePrecision.percentInCategory);
    numberOfPeople = 6;
    
   % drawPrecisionStatistics(meanPerPerson, stdPerPerson, randomSizePrecision.sizesOfClasses);
        
    experimentsSubjectMatrixLogical = false(size(experimentsSubjectMatrix,1), numberOfPeople);
    for i =1:numberOfPeople
        experimentsSubjectMatrixLogical(:,i) = experimentsSubjectMatrix == i;
    end
    
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
    
    % mean and std over people
    regionScores = mean(categoriesScores,2);
    regionStdScores = std(categoriesScores,1,2);
    
    regionColors = createColorMap(length(regionNames));
    figure(1234123);
    drawBars(regionScores,regionStdScores,regionNames,regionColors,'between / within');

    %find genes which appear both in GO and in Allen
    [commonGenes,commonGeneInGo,commonGeneInAllen] = intersect(geneGo.geneNames, selectedProbesData.gene_symbols);
    go_gene_mat = geneGo.go_genes_mat(:, commonGeneInGo);
    geneNames = geneGo.geneNames(commonGeneInGo);
    go_cat_without_genes = sum(go_gene_mat,2) == 0;
    cat_ids = geneGo.cat_ids(~go_cat_without_genes);
    go_gene_mat = go_gene_mat(~go_cat_without_genes, :);
    experimentsDataMatrix = experimentsDataMatrix(:, commonGeneInAllen);
    
    expressionDistances = distanceUsingMatrices(experimentsDataMatrix);
    
    diagIndices = 1:size(expressionDistances,1)+1:numel(expressionDistances); 
    expressionDistances(diagIndices) = inf;
    
    
    
    % for each category calc the new distance Score
    numOfGenesInCategory = full(sum(go_gene_mat,2));
    
    catDistanceScore = nan(length(cat_ids),numberOfPeople,numberOfRegions);
    
    
    
    for i = 1:size(catDistanceScore,1)
            currentCategorySubset = logical(full(go_gene_mat(i,:)));
            currentDistanceBetweenGroups = distanceBetweenGroups(:,:,:,~currentCategorySubset);
            %currentDistanceBetweenGroups = distanceBetweenGroups(:,:,:,currentCategorySubset);
            
            for m =1:numberOfRegions
                currentAreaDistanceBetweenGroups = squeeze(currentDistanceBetweenGroups(m,:,:,:));
                catDistanceScore(i,:,m) = clusteringScore( currentAreaDistanceBetweenGroups );
            end
            printPercentCounter(i, length(cat_ids) );
    end
    
%         
%     for m =1:numberOfRegions
%         fprintf('\n%s:   ', regionNames{m});
%         currentAreaSamplesIndices = experimentRegion(:,m);
%         currentAreaDataMatrix = experimentsDataMatrix(currentAreaSamplesIndices,:);
%         currentAreaExpressionDistances = distanceUsingMatrices(currentAreaDataMatrix);
%     
%         diagIndices = 1:size(currentAreaExpressionDistances,1)+1:numel(currentAreaExpressionDistances); 
%         currentAreaExpressionDistances(diagIndices) = inf;
%         currentAreaSamplesSubjectLabels = experimentsSubjectMatrixLogical(currentAreaSamplesIndices,:);
%         for i = 1:size(catDistanceScore,1)
%             currentCategorySubset = logical(full(go_gene_mat(i,:)));
%             newDistanceMatrix = removeFeaturesFromDistance(currentCategorySubset, currentAreaExpressionDistances, currentAreaDataMatrix);
%             catDistanceScore(i,:,m) = calcPrecisionAt1( newDistanceMatrix, currentAreaSamplesSubjectLabels, true);
%             
%             printPercentCounter(i, length(cat_ids) );
%         end
%     end
    
    save('catScores.mat', 'catDistanceScore', 'cat_ids', 'go_gene_mat','geneNames','experimentsDataMatrix','numOfGenesInCategory','regionNames','experimentRegion');
    
%    subsampleAndCheckPrecision(a, 'newRandomClasses.mat');
%    [a,b,c] = unique(numOfGenesInCategory);
end

function drawPrecisionStatistics(meanPerPerson, stdPerPerson, sizesOfClasses)

    h = figure(1111);
    plot(sizesOfClasses, stdPerPerson','.');
    xlabel('number of genes removed','fontsize',20)
    ylabel('std','fontsize',20);
    legend([repmat('human',[6,1]), num2str([1:6]')  ]); 
    h = figure(2222);
    plot(sizesOfClasses, meanPerPerson','.');
    xlabel('number of genes removed','fontsize',20)
    ylabel('mean','fontsize',20);
    legend([repmat('human',[6,1]), num2str([1:6]')  ],'fontsize',20); 
end

function [meanPerPerson, stdPerPerson] = calcNormalForEachPerson(percentInCategory)
    numberOfCategoriesSizes = size(percentInCategory,3);
    numberOfPeople = size(percentInCategory,2);
    meanPerPerson = nan(numberOfPeople, numberOfCategoriesSizes);
    stdPerPerson = nan(numberOfPeople, numberOfCategoriesSizes);
    for i =1: numberOfPeople
        for j =1: numberOfCategoriesSizes
            meanPerPerson(i,j) = mean(percentInCategory(:,i,j));
            stdPerPerson(i,j) = std(percentInCategory(:,i,j));
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

    distances = distanceUsingMatrices(onlyTheRemovedFeatures);
    
    newDistanceMatrix = distanceMatrix.^2 - distances.^2;

end

function distances = distancesUsingPdist(x)
        distanceVector = pdist(x,'euclidean');
        distances = squareform(distanceVector);
end

function distances = distanceUsingMatrices(x)
    Qx=repmat(dot(x,x,2),1,size(x,1));
    squaredDistances = (Qx+Qx'-2*(x*x'));
    distances = sqrt(squaredDistances);
end
