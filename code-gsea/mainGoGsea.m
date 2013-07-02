
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

goCategoriesGseaScores = cell(size(anovaResults.regionLabels));
goCategoriesGseaPvalScores = cell(size(anovaResults.regionLabels));

s = RandStream('mt19937ar','Seed',0);
no_of_permutations=1000;
randPermMatrix = nan(length(geneGoNames), no_of_permutations);
for i_perm = 1:no_of_permutations
   randPermMatrix(:,i_perm) = randperm(s,length(geneGoNames));
%   randPermMatrix(:,i_perm) = randperm(length(geneGoNames));
end
%tic

%replace zero with the smallest number
zerosEntries = anovaMatrix(:) == 0;
anovaMatrix(zerosEntries) = min(anovaMatrix(~zerosEntries));

absoulteLogOfAnovaScroes = -log(anovaMatrix);

for i =1:length(anovaResults.regionLabels)
   fprintf('\n %s:   ', anovaResults.regionLabels{i});
   regionAnovaScores = absoulteLogOfAnovaScroes(:,i);
   [sortedGeneByAnovaScore, sortIndices] = sort( regionAnovaScores);
   sortedCat2GeneTable = go_genes_mat(:,sortIndices);
   
   goCatScroes = nan(size(cat_ids));
   goCatPvalScores = nan(size(cat_ids));
%    tic
   for j=1:length(cat_ids)
       currentGoCatMembers = sortedCat2GeneTable(j,:)';
       [pval,~] = GSEA(sortedGeneByAnovaScore,currentGoCatMembers,1,randPermMatrix);
       %goCatScroes(j) = es;
       goCatPvalScores(j) = pval;
       printPercentCounter(j, length(cat_ids));
%        toc
   end
   %goCategoriesGseaScores{i} = goCatScroes;
   goCategoriesGseaPvalScores{i} = goCatPvalScores;
end
%save('goCatScores.mat');
analyzeCatWithScores(cat_ids, aspects, goCategoriesGseaPvalScores, regionLabels,'houseKeeping.html','rank sum results',{'house keeping genes'} );