function [mapping ,genesNotFound] = createAllenToEntrezMapping(listOfGenes,dataType)
% mapping = createAllenToEntrezMapping(listOfGenes,dataType)
% maps allen_gene_names <-> entrez <-> ncbi_access_number
% possible datatypes: 'allen' , 'entrez' 
% 
% mapping = [entrez, ncbi_access_number, allen_gene_symbol, allen_long_gene_name]


%[geneid	genesymbol	genename	imageseriesid	plane	entrezgeneid	ncbiaccessionnumber] = textread('gene-series-map.tsv', '%s %s %s %s %s %s %s','delimiter','\t','headerlines',1);

% [~,	genesymbol,	genename,	~,	~,	entrezgeneid,	ncbiaccessionnumber] = textread('gene-series-map.tsv', '%s %s %s %s %s %d %s','delimiter','\t','headerlines',1);
%     entrez = entrezgeneid;
%     entrez( entrezgeneid == 0) = nan;
% save('Allen2Entrez.mat');

switch dataType
    case 'allen'
        [mapping ,genesNotFound] = allen2Entrez(listOfGenes);
    case 'entrez'
        [mapping ,genesNotFound] = entrez2Allen(listOfGenes);
end


end


function [mapping ,genesNotFound]= allen2Entrez(gene_names)

    load('Allen2Entrez.mat');
%     [foundInList, geneInd] = ismember(gene_names, genesymbol);
    geneInd = cellfind(genesymbol, gene_names);
    genesNotFound = geneInd == 0;
    releventInd = geneInd(~genesNotFound);
    
    mapping = cell(length(gene_names), 4);
    mapping(~genesNotFound,:) = [num2cell(entrez(releventInd)), ncbiaccessionnumber(releventInd), genesymbol(releventInd) ,genename(releventInd)];
end

function entrez2Allen(entrez_input)

    load('Allen2Entrez.mat');
    
   
end