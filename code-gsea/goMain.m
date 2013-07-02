% GO = geneont('live',true);
%GO = geneont('file', 'gene_ontology.obo.txt');
filename = 'gene_association.goa_human';
genes = {'Calb1' , 'DDX3Y','TTTY15'};

[go_genes_mat, cat_ids, geneNames, aspects] = get_gos_for_all_genes(filename, GO);

%save('humanGene2GoMatrix.mat', 'go_genes_mat', 'cat_ids', 'geneNames', 'aspects');
[go_genes_mat, cat_ids, aspects] = get_gos_for_genes(genes, filename, GO);


load('humanGene2GoMatrix.mat', 'go_genes_mat', 'cat_ids', 'geneNames', 'aspects');

go_cat_indices = cat_ids(logical(go_genes_mat(:,1)));
names = get(GO(go_cat_indices).terms, 'name');

%go_cat_names = go_id2name(go_cat_indices, GO);

[pval,ES] = inHouse_GSEA(sorted_pred,sorted_labels,pwr);