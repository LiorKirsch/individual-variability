function [go_genes_mat, cat_ids, geneNames, aspects] = get_gos_for_all_genes(filename, GO)
% this function outputs a matrix, num_go_categories* num_genes, which shows
% for each go category the genes of interest associated with
% it.
% inputs:
% genes - cell array of gene names (strings)
% gos - vector of go categories of interest (double)
% filename - filename of GO for animal of interest
% GO - updated GO structure
% 
% HGann = goannotread(filename);
% 
% HGmap = containers.Map();
% for i=1:numel(HGann)
%     key = HGann(i).DB_Object_Symbol;
%     if isKey(HGmap,key)
%         HGmap(key) = [HGmap(key) HGann(i).GOid];
%     else
%         HGmap(key) = HGann(i).GOid;
%     end
% end
% 
% end

% this function outputs a matrix, num_go_categories* num_genes, which shows
% for each go category the genes of interest associated with
% it.
% inputs:
% genes - cell array of gene names (strings)
% gos - vector of go categories of interest (double)
% filename - filename of GO for animal of interest
% GO - updated GO structure


% get GO annotations for animal of interest
%SGDann = goannotread_myfix(filename);
SGDann = goannotread(filename);

SGDgenes  = {SGDann.DB_Object_Symbol}; %  gene list
SGDgo     = [SGDann.GOid];             % associated GO terms
SGDaspect = {SGDann.Aspect};

[all_terms, idx1, ~] = unique(SGDgo);
aspects = SGDaspect(idx1);
num_terms = length(all_terms);

[uniqueGenes, ~, geneIinGo] = unique(SGDgenes);
num_genes = length(uniqueGenes);

%go_cats = cell(num_genes,1);
%go_count = zeros(num_terms,1);    % a vector of GO term counts for the entire chip.
go_genes_mat = (sparse(num_terms,num_genes));

fprintf('Creating genes to Go metrix:    ');
for i_gene = 1:num_genes
    %idx = strcmpi(SGDgenes,uniqueGenes{i_gene});  % lookup gene
    goid = SGDgo(geneIinGo == i_gene);                 % get the respective GO ids

    relevant_terms = getancestors(GO,goid);
    
    %get The Index Of The ReleveantTerms in the all_terms
    [~, term_inds] = intersect(all_terms, relevant_terms);
    
    printPercentCounter(i_gene, num_genes);
    go_genes_mat(term_inds, i_gene) = 1;
end

cat_ids = all_terms;
geneNames = uniqueGenes;
go_genes_mat = logical(go_genes_mat);
end