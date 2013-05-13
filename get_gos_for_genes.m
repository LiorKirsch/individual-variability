function [go_genes_mat, cat_ids, aspects] = get_gos_for_genes(genes, filename, GO)
% this function outputs a matrix, num_go_categories* num_genes, which shows
% for each go category the genes of interest associated with
% it.
% inputs:
% genes - cell array of gene names (strings)
% gos - vector of go categories of interest (double)
% filename - filename of GO for animal of interest
% GO - updated GO structure

num_genes = length(genes);

% get GO annotations for animal of interest
%SGDann = goannotread_myfix(filename);
SGDann = goannotread(filename);

%SGDann = goannotread(filename);
SGDgenes  = {SGDann.DB_Object_Symbol}; %  gene list
SGDgo     = [SGDann.GOid];             % associated GO terms
SGDaspect = {SGDann.Aspect};

[all_terms, idx1, junk] = unique(SGDgo);
aspects = SGDaspect(idx1);
num_terms = length(all_terms);

%go_cats = cell(num_genes,1);
%go_count = zeros(num_terms,1);    % a vector of GO term counts for the entire chip.
go_genes_mat = (zeros(num_terms,num_genes));
t1 = clock;
for i_gene = 1:num_genes
    t2 = clock;
    fprintf('Iteration %d out of %d, time elapsed %g sec \n', i_gene, num_genes, etime(t2, t1))
    idx = strcmpi(SGDgenes,genes{i_gene});  % lookup gene
    goid = SGDgo(idx);                 % get the respective GO ids

    relevant_terms = getancestors(GO,goid);
    [junk, term_inds] = intersect(all_terms, relevant_terms);

    go_genes_mat(term_inds, i_gene) = 1;
end

cat_ids = all_terms;
go_genes_mat = sparse(go_genes_mat);