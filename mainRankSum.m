anovaResults = load('anovaResults.mat','anovaMatrix', 'selectedProbesData', 'regionLabels','samplesInRegion');

%geneGo = load('humanGene2GoMatrix.mat', 'go_genes_mat', 'cat_ids', 'geneNames', 'aspects');
geneGo = load('housekeepingGenes.mat', 'go_genes_mat', 'cat_ids', 'geneNames', 'aspects');

probesGeneNames = anovaResults.selectedProbesData.gene_symbols;

[jointGenes, goIntersectIndices,probeIntersectIndices] = intersect(geneGo.geneNames, probesGeneNames );

probesGeneNames = probesGeneNames(probeIntersectIndices);
anovaMatrix = anovaResults.anovaMatrix(probeIntersectIndices,:);
geneGoNames = geneGo.geneNames(goIntersectIndices);
go_genes_mat = geneGo.go_genes_mat(:,goIntersectIndices);

% remove regions that had less then 10 smaples (mean across subjects)
areaThreshold = mean(anovaResults.samplesInRegion,2) >= 10;
anovaMatrix = anovaMatrix(:, areaThreshold);
regionLabels = anovaResults.regionLabels(areaThreshold);
clear('areaThreshold');

% also remove categories with no genes
numberOfGenesInCategory = full(sum(go_genes_mat,2));
go_genes_mat = go_genes_mat(~(numberOfGenesInCategory==0),:);
go_genes_mat = logical(go_genes_mat);
cat_ids = geneGo.cat_ids(~(numberOfGenesInCategory==0));
aspects = geneGo.aspects(~(numberOfGenesInCategory==0));
    
goCategoriesGseaScores = cell(size(regionLabels));
goCategoriesGseaPvalScores = cell(size(regionLabels));


%replace zero with the smallest number
zerosEntries = anovaMatrix(:) == 0;
anovaMatrix(zerosEntries) = min(anovaMatrix(~zerosEntries));

Log10AnovaScores = -log10(anovaMatrix);

for i_region =1:length(regionLabels)
   fprintf('\n %s:   ', regionLabels{i_region});
   regionAnovaScores = Log10AnovaScores(:, i_region);
   [sortedGeneByAnovaScore, sortIndices] = sort(regionAnovaScores);
   sortedCat2GeneTable = go_genes_mat(:,sortIndices);
   
   goCatPvalScores = nan(size(cat_ids));

   for j=1:length(cat_ids)
       smallGroup =  sortedCat2GeneTable(j,:)';
       largeGroup = ~smallGroup;
       [p, h, stats] = ranksum(sortedGeneByAnovaScore(smallGroup), sortedGeneByAnovaScore(largeGroup) ); 
       goCatPvalScores(j) = p;
       printPercentCounter(j, length(cat_ids));
   end
   %goCategoriesGseaScores{i_region} = goCatScroes;
   %fdrCorrectedValues = mafdr(goCatPvalScores, 'BHFDR', true);
   goCategoriesGseaPvalScores{i_region} = goCatPvalScores;
end
%save('ranksumScores.mat');
%analyzeCatWithScores(cat_ids, aspects, goCategoriesGseaPvalScores, regionLabels,'ranksum.html','rank sum results');
analyzeCatWithScores(cat_ids, aspects, goCategoriesGseaPvalScores, regionLabels,'houseKeeping.html','rank sum results',{'house keeping genes'} );