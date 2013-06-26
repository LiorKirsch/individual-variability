function subWithGoAndCheckPrecision()

    close all;

    geneGo = load('humanGene2GoMatrix.mat', 'go_genes_mat', 'cat_ids', 'geneNames', 'aspects');
    
    %load('dataWithDistance.mat','experimentsSubjectMatrix', 'experimentsDataMatrix', 'expressionDistances','experimentRegion','regionNames');
    %load('expressionRawDataWithOntology.mat', 'selectedProbesData');
    
    load('pcaData.mat','experimentsSubjectMatrix', 'experimentsDataMatrix', 'experimentRegion','regionNames','selectedProbesData');
%     randomSizePrecision =  load('randomSizePrecision.mat','percentInCategory','sizesOfClasses');
%     
%     [meanPerPerson, stdPerPerson] = calcNormalForEachPerson(randomSizePrecision.percentInCategory);
    numberOfPeople = size(experimentsSubjectMatrix,2);
    
   % drawPrecisionStatistics(meanPerPerson, stdPerPerson, randomSizePrecision.sizesOfClasses);
        
    numberOfSamples = size(experimentsSubjectMatrix,1);
    numberOfRegions = size(experimentRegion,2);
    assert(  numberOfSamples == size(experimentsDataMatrix,1) );

    
    [categoriesScores, ~] = calcRegionBetweenSubjectDist(experimentsSubjectMatrix, experimentsDataMatrix, experimentRegion);
    
    % mean and std over people
    regionScores = mean(categoriesScores,2);
    regionStdScores = std(categoriesScores,1,2);
    regionColors = createColorMap(length(regionNames));
    figure(1234123);
    drawBars(regionScores,regionStdScores,regionNames,regionColors,'between / within');

%     addpath('~/Projects/general use functions/classifiy to neuron astro oligo/');
%     entrez_ids = selectedProbesData.entrez_ids;
%     neuronGliaClass = cahoyClassifyGene(entrez_ids, true);
%     go_gene_mat = neuronGliaClass';
%     cat_ids = {'astrocytes','neurons','oligodendrocytes'};
    
    %find genes which appear both in GO and in Allen
    [commonGenes,commonGeneInGo,commonGeneInAllen] = intersect(geneGo.geneNames, selectedProbesData.gene_symbols);
    entrez_ids = selectedProbesData.entrez_ids( commonGeneInAllen );
    experimentsDataMatrix = experimentsDataMatrix(:, commonGeneInAllen);
    go_gene_mat = geneGo.go_genes_mat(:, commonGeneInGo);
    geneNames = geneGo.geneNames(commonGeneInGo);
    
    selectedProbesData.probe_ids = selectedProbesData.probe_ids(commonGeneInAllen);
    selectedProbesData.probe_names = selectedProbesData.probe_names(commonGeneInAllen);
    selectedProbesData.gene_ids = selectedProbesData.gene_ids(commonGeneInAllen);
    selectedProbesData.gene_symbols = selectedProbesData.gene_symbols(commonGeneInAllen);
    selectedProbesData.gene_names = selectedProbesData.gene_names(commonGeneInAllen);
    selectedProbesData.entrez_ids = selectedProbesData.entrez_ids(commonGeneInAllen);
    selectedProbesData.chromosome = selectedProbesData.chromosome(commonGeneInAllen);
    selectedProbesData.bestProbeForGene = selectedProbesData.bestProbeForGene(commonGeneInAllen);
    
    go_cat_without_genes = sum(go_gene_mat,2) == 0;
    cat_ids = geneGo.cat_ids(~go_cat_without_genes);
    go_gene_mat = go_gene_mat(~go_cat_without_genes, :);
    aspects = geneGo.aspects(~go_cat_without_genes);
    
%     houseKeepingFileName = '/home/lab/ossnat/work/Data/Homo/Expression/HousekeepingGenes/Erez_2003_RefSeq_mRNA.txt';
%     [~, ~, houseKeepingEntrezGeneID, ~] = textread(houseKeepingFileName, '%s %s %d %s', 'headerlines', 1,'delimiter','\t');
%     go_gene_mat = ismember(entrez_ids, houseKeepingEntrezGeneID)';
%     cat_ids = {'house keeping genes'};
%      
    % for each category calc the new distance Score
    numOfGenesInCategory = full(sum(go_gene_mat,2));
    
    
    [categoriesScores, distanceBetweenGroups] = calcRegionBetweenSubjectDist(experimentsSubjectMatrix, experimentsDataMatrix, experimentRegion);
        % mean and std over people
    regionScores = mean(categoriesScores,2);
    regionStdScores = std(categoriesScores,1,2);
    regionColors = createColorMap(length(regionNames));
    figure(12234123);
    drawBars(regionScores,regionStdScores,regionNames,regionColors,'between / within');

    
    catDistanceScore = nan(length(cat_ids),numberOfPeople,numberOfRegions);
    onlyCatDistanceScore = nan(length(cat_ids),numberOfPeople,numberOfRegions);

    
    fprintf('\nrunning over categories:   ');
    for i = 1:size(catDistanceScore,1)
            currentCategorySubset = logical(full(go_gene_mat(i,:)));
            currentDistanceBetweenGroups = distanceBetweenGroups(:,:,:,~currentCategorySubset);
            currentDistanceBetweenGroupsOnlyCat = distanceBetweenGroups(:,:,:,currentCategorySubset);
            
            parfor m =1:numberOfRegions
                currentAreaDistanceBetweenGroups = squeeze(currentDistanceBetweenGroups(m,:,:,:));
                catDistanceScore(i,:,m) = clusteringScore( currentAreaDistanceBetweenGroups );
                
                currentAreaOnlyCatDistanceBetweenGroups = squeeze(currentDistanceBetweenGroupsOnlyCat(m,:,:,:));
                onlyCatDistanceScore(i,:,m) = clusteringScore( currentAreaOnlyCatDistanceBetweenGroups );
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
    
     save('catScores.mat', 'catDistanceScore', 'onlyCatDistanceScore','cat_ids', 'aspects','go_gene_mat','geneNames','experimentsDataMatrix','numOfGenesInCategory','regionNames','experimentRegion','selectedProbesData');
%     save('catHouseKeepingScores.mat', 'catDistanceScore', 'onlyCatDistanceScore','cat_ids', 'go_gene_mat','geneNames','experimentsDataMatrix','numOfGenesInCategory','regionNames','experimentRegion');
    
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
