 serotonin_pathway_genes = { 'HTR2A','HTR2B','HTR2C',...
        'HTR1A','HTR1B','HTR1D','HTR1E','HTR1F',...
            'HTR3A','HTR3B','HTR3C','HTR3D','HTR3E',...
            'HTR4','HTR5A','HTR6','HTR7',...
            'TPH2','TPH1','DDC',...
            'SLC18A1','SLC18A2',...
            'SLC6A4','MAOA','MAOB'};

dopamin_pathway_genes = {'DRD1','DRD2','DRD3','DRD4','DRD5','SLC6A3','MAOA','MAOB','COMT',...
            'TH','SLC18A1','SLC18A2','DDC'};
        
        

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
go_genes_mat = logical(go_genes_mat);
cat_ids = geneGo.cat_ids(~(numberOfGenesInCategory==0));

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


[~,set_indices] = ismember(serotonin_pathway_genes, geneGoNames);
seratoninBolleanList = false(size(geneGoNames))';
seratoninBolleanList(set_indices) = true;

[~,set_indices] = ismember(dopamin_pathway_genes, geneGoNames);
dopaminBolleanList = false(size(geneGoNames))';
dopaminBolleanList(set_indices) = true;

sertGseaPvalScores = nan(length(anovaResults.regionLabels),1);
dopaminGseaPvalScores = nan(length(anovaResults.regionLabels),1);

for i =1:length(anovaResults.regionLabels)
   fprintf('\n %s:   ', anovaResults.regionLabels{i});
   regionAnovaScores = absoulteLogOfAnovaScroes(:,i);
   [sortedGeneByAnovaScore, sortIndices] = sort( regionAnovaScores);
   sortedCat2GeneTable = go_genes_mat(:,sortIndices);
   
   [pval,es] = GSEA(sortedGeneByAnovaScore,seratoninBolleanList,1,randPermMatrix);
   sertGseaPvalScores(i) = pval;
   
   [pval,~] = GSEA(sortedGeneByAnovaScore,dopaminBolleanList,1,randPermMatrix);
   dopaminGseaPvalScores(i) = pval;
end

scores = [dopaminGseaPvalScores , sertGseaPvalScores]';

catLabels = repmat({{'seratonin pathway', 'dopamin pathway'}}, [length(anovaResults.regionLabels),1]);
printCatScoreIntoHtml(anovaResults.regionLabels, scores, scores, catLabels, 'seratonin and dopamin pathways', 'sertDop.html' )

save('sertDopScores.mat');
