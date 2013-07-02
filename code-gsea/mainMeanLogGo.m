
anovaResults = load('anovaResults.mat','anovaMatrix', 'selectedProbesData', 'regionLabels');


geneGo = load('humanGene2GoMatrix.mat', 'go_genes_mat', 'cat_ids', 'geneNames', 'aspects');

probesGeneNames = anovaResults.selectedProbesData.gene_symbols;

[jointGenes, goIntersectIndices,probeIntersectIndices] = intersect(geneGo.geneNames, probesGeneNames );

probesGeneNames = probesGeneNames(probeIntersectIndices);
anovaMatrix = anovaResults.anovaMatrix(probeIntersectIndices,:);
geneGoNames = geneGo.geneNames(goIntersectIndices);
go_genes_mat = geneGo.go_genes_mat(:,goIntersectIndices);

% also remove categories with no genes
numberOfGenesInCategory = full(sum(go_genes_mat,2));
go_genes_mat = go_genes_mat(~(numberOfGenesInCategory==0),:);
cat_ids = geneGo.cat_ids(~(numberOfGenesInCategory==0));

goCategoriesMeanLogScores = nan(length(cat_ids), length(anovaResults.regionLabels));


%replace zero with the smallest number
zerosEntries = anovaMatrix(:) == 0;
anovaMatrix(zerosEntries) = min(anovaMatrix(~zerosEntries));
logAnovaMatrix = log(anovaMatrix);
numberOfGenesInCategory = full(sum(go_genes_mat,2));


for i =1:length(anovaResults.regionLabels)
   regionAnovaScores = logAnovaMatrix(:,i);
   
   goCatMembersScoresSum = go_genes_mat *regionAnovaScores;
   meanLogPScore = goCatMembersScoresSum ./ numberOfGenesInCategory;
%    for j=1:length(geneGo.cat_ids)
%        currentGoCatMembers = sortedCat2GeneTable(j,:);
%        goCatMembersScroes = regionAnovaScores(logical(currentGoCatMembers) );
%        meanLogPScore = mean(log( goCatMembersScroes ) );
%        
%        goCatPvalScores(j) = meanLogPScore;
%        
%    end
   %goCategoriesGseaScores{i} = goCatScroes;
   backToPvalues = exp(meanLogPScore);
 %  fdrCorrectedValues = mafdr(backToPvalues, 'BHFDR', true);
   [pthr,pcor,padj]  = fdr(backToPvalues);

   goCategoriesMeanLogScores(:,i) = pcor;
end

sortedValues = nan(size(goCategoriesMeanLogScores));
sortedCatIndices = nan(size(goCategoriesMeanLogScores));
sortedCatNames = cell(size(goCategoriesMeanLogScores,2),1);

 GO = geneont('file', 'gene_ontology.obo.txt');
 
catNames = go_id2name(cat_ids, GO);
for i =1:length(anovaResults.regionLabels)
   [sortedValues(:,i), sortIndices] = sort(goCategoriesMeanLogScores(:,i));
   sortedCatIndices(:,i) = cat_ids(sortIndices);
   sortedCatNames{i} = catNames(sortIndices);

end

printCatScoreIntoHtml(anovaResults.regionLabels, sortedValues, sortedCatIndices, sortedCatNames, 'mean log scores', 'meanLog.html' )
