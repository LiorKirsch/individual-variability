function [cats, pvals, out_aspects, num_genes_in_clust, cat_gene_names] = ...
    GO_enrich_ABA(clust_genes, all_genes, GO, db, tree, brain_related, mc, plot)

% this part loads a pre-made matrix of genes * go categories for several
% datasets that I use often. every once in a while, I'll update the
% matrix. I produce these matrices using the function get_gos_for_genes.
if strcmp(db, 'adult')
    load('go_genes_mat_adult_mouse')
    gene_names = aba_gene_names;
elseif strcmp(db, 'mouse')  % this is the relevant matrix for the mouse-hourglass project.
    %it was created only for the 2002 genes in the dataset, and here it is
    %screened for the background genes (the input all_genes) only.
    load('dev_go_genes_mat')
    ginds = cellfind(gene_names, all_genes);
    go_genes_mat = go_genes_mat(:, ginds);
    gene_names = gene_names(ginds);
elseif strcmp(db, 'human_affy')
    load('h_affy_go_genes_mat')
elseif strcmp(db, 'human_rseq')
    load('h_rseq_go_genes_mat')
elseif strcmp(db, 'cb')
    load('cb_go_genes_mat')
    ginds = cellfind(gene_names, all_genes);
    go_genes_mat = go_genes_mat(:, ginds);
    gene_names = gene_names(ginds);
else
    error('unknown database')
end

if strcmp(brain_related, 'yes')
%     [go_genes_mat, cat_ids, aspects, cat_names] = ...
%         screen_non_brain_related_cats(go_genes_mat, cat_ids, aspects, cat_names);
% 
load('images_go_genes_mat_screened', 'brain_cat_ids')
[~, relevant_inds, ~] = intersect(cat_ids, brain_cat_ids);
go_genes_mat = go_genes_mat(relevant_inds, :);
cat_names = cat_names(relevant_inds);
cat_ids = cat_ids(relevant_inds);
aspects = aspects(relevant_inds);
end

% this is where I set the significance cutoff
alpha = 0.1;

% screening the go categories for those of a specific DAG (molecular
% function, biological process or cellular component, if specified by the
% input parameter 'tree'.
if ~(strcmp(tree, 'all'))
    inds = cellfind(aspects, tree);
    go_genes_mat(setdiff([1:size(go_genes_mat, 1)], inds), :) = [];
    cat_ids = cat_ids(inds);
    aspects = aspects(inds);
end


% screen for GO categories with number of genes above a certain cutoff
cutoff = 10;
num_genes_in_cat = sum(go_genes_mat, 2);
relevant_inds = find(num_genes_in_cat>=cutoff);
go_genes_mat = go_genes_mat(relevant_inds, :);
cat_ids = cat_ids(relevant_inds);
aspects = aspects(relevant_inds);

% here is where we search for enrichment
all_terms = cat_ids;
num_terms = length(all_terms);

all_inds = cellfind(gene_names, all_genes);
all_inds = all_inds(find(all_inds));
clust_inds = cellfind(gene_names, clust_genes);
clust_inds = clust_inds(find(clust_inds));

all_genes_count = sum(go_genes_mat(:, all_inds), 2);
clust_genes_count = sum(go_genes_mat(:, clust_inds), 2);

%pvalues = hygepdf(clust_genes_count, sum(all_genes_count), all_genes_count, sum(clust_genes_count));
pvalues = hygepdf(clust_genes_count, length(all_inds), all_genes_count, length(clust_inds));

%this is for determining whether the category is under represented or over
%reresented:
dir = sign(clust_genes_count/sum(clust_genes_count) - all_genes_count/sum(all_genes_count));
dir = sign(clust_genes_count/length(clust_inds) - all_genes_count/length(all_inds));

% correct for multiple comparisons:
switch mc
    case 'bonf'
        num_comparisons = num_terms;
        cutoff = alpha / num_comparisons;
        inds1 = find(pvalues <= cutoff);
        inds2 = find(dir==1);
        inds = intersect(inds1, inds2);
        pvals(:, 1) = pvalues(inds);
        pvals(:,2) = pvalues(inds)*num_comparisons;
    case 'fdr'
        num_comparisons = num_terms;
        [pvaluesBH] = mafdr(pvalues, 'BHFDR', true);
        %[pvaluesBH] = mafdr(pvalues);

        inds1 = find(pvaluesBH <= alpha);
        inds2 = find(dir==1);
        inds = intersect(inds1, inds2);
        pvals(:, 1) = pvalues(inds);
        pvals(:,2) = pvaluesBH(inds);
        pvals(:,3) = pvalues(inds)*num_comparisons;
    case 'none'
        inds1 = find(pvalues <= alpha);
        inds2 = find(dir==1);
        inds = intersect(inds1, inds2);
        pvals(:, 1) = pvalues(inds);
    otherwise
        error('wrong type of multiple comparison correction')
end

% getting the output parameters: the category ids, their aspects (useful if
% we didn't specify a certain DAG), the number of genes in the cluster and
% the gene names of that category in the cluster.
cats = all_terms(inds);
out_aspects = aspects(inds);
num_genes_in_clust = clust_genes_count((inds));

cat_gene_names = cell(length(cats), 1);
for i_cat = 1:length(cats)
    cat_gene_names{i_cat} = intersect(gene_names(find(go_genes_mat((cat_ids==cats(i_cat)), :))), clust_genes);
end

% some output - plot an enrichment tree if needed. I havent used this in a
% long long time and not sure if it works. also this part outputs whether
% there's enrichment found.
if ~isempty(pvals)
    disp('Enrichment found!')
    if strcmp(plot, 'plot')
        % visualize categories in a tree
        subGO = GO(getancestors(GO,cats));
        [cm acc rels] = getmatrix(subGO);
        BG = biograph(cm,get(subGO.Terms,'name'));
        % Use the p-values, calculated before, to assign a color to the graph
        % nodes. In this example an arbitrary color map is used, where bright red
        % is the most significant and bright green is the least significant.
        num_genes = zeros(length(acc),1);
        for i=1:numel(acc)
            pval = pvalues(find(cat_ids==acc(i)));
            if ~isempty(pval)
                color = [(1-pval).^(10),pval.^(1/10),0.3];
                set(BG.Nodes(i),'Color',color);
            end
            current_cat = go_id2name(acc(i), GO);
            num_genes(i) = clust_genes_count(cat_ids==acc(i));
            set(BG.Nodes(i),'Label',sprintf('%s\n%3.1d\n#genes = %d',current_cat{:}, pval, num_genes(i))) % add info to datatips
        end
        %view(BG);
        num_genes_in_clust = num_genes;
        BG.view;
        % save a figure of the tree
        %print_bg_to_figure(BG)
    end
else
    disp('No enriched categories')
end
end
%end