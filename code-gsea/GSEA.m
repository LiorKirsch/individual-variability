function [pval,max_ES_1] = GSEA(sorted_pred,sorted_labels,pwr,randPermMatrix)
% lior kirsch 2013 based on Ossnat mar. 2011

assert(all(size(sorted_pred) == size(sorted_labels)));
%assert( islogical(sorted_labels) );

%no_of_permutations=1000;
no_of_permutations=size(randPermMatrix,2);

[max_ES_1,ES] = calc_max_ES2(sorted_pred,sorted_labels,pwr);

%addpath('/private/biu/packages/matlab/R2009b/toolbox/matlab/randfun');

%  permuted_labels_matrix = sorted_labels(randPermMatrix);
EM_perms = nan(no_of_permutations,1);
% create permutation of labels
parfor i_perm = 1:no_of_permutations
    permVector = randPermMatrix(:,i_perm);
    permuted_labels = sorted_labels(permVector); 
    EM_perms(i_perm) = calc_max_ES2(sorted_pred,permuted_labels,pwr);
    
%     EM_perms(i_perm) = calc_max_ES2(sorted_pred,permuted_labels_matrix(:,i_perm),pwr);
end


% hist(EM_perms,no_of_permutations)
% DBG parms - 
mean_ = mean(EM_perms);
std_ = std(EM_perms);
max_ = max(EM_perms);

% p-value - calc Z score

%addpath('/usr/local/matlab2008a/toolbox/matlab/specfun')
[h,pval,ci,zval] = ztest(max_ES_1,mean(EM_perms),std(EM_perms));


end
