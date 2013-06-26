function clacScoreWithoutAsingleGene()

    close all;

    geneGo = load('humanGene2GoMatrix.mat', 'go_genes_mat', 'cat_ids', 'geneNames', 'aspects');
    
    load('pcaData.mat','experimentsSubjectMatrix', 'experimentsDataMatrix', 'expressionDistances','experimentRegion','regionNames','selectedProbesData');
    
    numberOfPeople = size(experimentsSubjectMatrix,2);
    numberOfSamples = size(experimentsSubjectMatrix,1);
    numberOfRegions = size(experimentRegion,2);
    assert(  numberOfSamples == size(experimentsDataMatrix,1) );

    %find genes which appear both in GO and in Allen
    [commonGenes,commonGeneInGo,commonGeneInAllen] = intersect(geneGo.geneNames, selectedProbesData.gene_symbols);
    entrez_ids = selectedProbesData.entrez_ids( commonGeneInAllen );
    experimentsDataMatrix = experimentsDataMatrix(:, commonGeneInAllen);
    go_gene_mat = geneGo.go_genes_mat(:, commonGeneInGo);
    geneNames = geneGo.geneNames(commonGeneInGo);
    
    go_cat_without_genes = sum(go_gene_mat,2) == 0;
    cat_ids = geneGo.cat_ids(~go_cat_without_genes);
    go_gene_mat = go_gene_mat(~go_cat_without_genes, :);
    aspects = geneGo.aspects(~go_cat_without_genes);
    
    numOfGenesInCategory = full(sum(go_gene_mat,2));
    numberOfGenes = length(geneNames);
    
    [categoriesScores, distanceBetweenGroups] = calcRegionBetweenSubjectDist(experimentsSubjectMatrix, experimentsDataMatrix, experimentRegion);
        % mean and std over people
      
    singleGeneDistanceScore = nan(numberOfGenes,numberOfPeople,numberOfRegions);
    onlySingleGeneDistanceScore = nan(numberOfGenes,numberOfPeople,numberOfRegions);

    
    fprintf('\nrunning over categories:   ');
    for i = 1:numberOfGenes
            currentCategorySubset = false(numberOfGenes,1);
            currentCategorySubset(i) = true;
            currentDistanceBetweenGroups = distanceBetweenGroups(:,:,:,~currentCategorySubset);
            currentDistanceBetweenGroupsOnlyCat = distanceBetweenGroups(:,:,:,currentCategorySubset);
            
            parfor m =1:numberOfRegions
                currentAreaDistanceBetweenGroups = squeeze(currentDistanceBetweenGroups(m,:,:,:));
                singleGeneDistanceScore(i,:,m) = clusteringScore( currentAreaDistanceBetweenGroups );
                
                currentAreaOnlyCatDistanceBetweenGroups = squeeze(currentDistanceBetweenGroupsOnlyCat(m,:,:,:));
                onlySingleGeneDistanceScore(i,:,m) = clusteringScore( currentAreaOnlyCatDistanceBetweenGroups );
            end
            printPercentCounter(i, numberOfGenes);
    end
  
     save('singleGeneScores.mat', 'singleGeneDistanceScore', 'onlySingleGeneDistanceScore','cat_ids', 'aspects','go_gene_mat','geneNames','experimentsDataMatrix','numOfGenesInCategory','regionNames','experimentRegion');
%     save('catHouseKeepingScores.mat', 'catDistanceScore', 'onlyCatDistanceScore','cat_ids', 'go_gene_mat','geneNames','experimentsDataMatrix','numOfGenesInCategory','regionNames','experimentRegion');
    
end

function distances = distanceUsingMatrices(x)
    Qx=repmat(dot(x,x,2),1,size(x,1));
    squaredDistances = (Qx+Qx'-2*(x*x'));
    distances = sqrt(squaredDistances);
end
