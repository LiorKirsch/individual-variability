function createHousekeepingList()

fileName = '/home/lab/ossnat/work/Data/Homo/Expression/HousekeepingGenes/Erez_2003_RefSeq_mRNA.txt';

[ensemblGeneID, ensemblTranscriptID, entrezGeneID, HGNCsymbol] = textread(fileName, '%s %s %d %s', 'headerlines', 1,'delimiter','\t');
clear('fileName');
% Ensembl Gene ID Ensembl Transcript ID   EntrezGene ID   HGNC symbol

anovaResults = load('anovaResults.mat','anovaMatrix', 'selectedProbesData', 'regionLabels','samplesInRegion');
anovaEntrez = anovaResults.selectedProbesData.entrez_ids;
probesGeneNames = anovaResults.selectedProbesData.gene_symbols;

[jointGenes, probeIntersectIndices,HGNCsymbolIndices] = intersect(anovaEntrez, entrezGeneID );
%[jointGenes2, probeIntersectIndices, HGNCsymbolIndices] = intersect(probesGeneNames, HGNCsymbol );



go_genes_mat = false(1,length(probesGeneNames));
go_genes_mat(probeIntersectIndices) = true;

cat_ids = [1];
geneNames = probesGeneNames;
aspects = {''};
save('housekeepingGenes.mat');

end