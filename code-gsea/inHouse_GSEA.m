function [pval,ES] = ...
    inHouse_GSEA(sorted_pred,sorted_labels,pwr)

% Ossnat mar. 2011

no_of_permutations=1000;
if min(sorted_labels)<0
    sorted_labels(find(sorted_labels<0))=0; % D B G - check met!!
end
if size(sorted_labels,1)~=size(sorted_pred,1)
    sorted_labels =sorted_labels';
end
[max_ES_1,ES] = calc_max_ES(sorted_pred,sorted_labels,pwr);

%addpath('/private/biu/packages/matlab/R2009b/toolbox/matlab/randfun');

% create permutation of labels
for i_perm = 1:no_of_permutations
    permuted_labels = sorted_labels(randperm(length(sorted_labels))); 
    EM_perms(i_perm) = calc_max_ES(sorted_pred,permuted_labels,pwr);
end


% hist(EM_perms,no_of_permutations)
% DBG parms - 
mean_ = mean(EM_perms)
std_ = std(EM_perms)
max_ = max(EM_perms)

% p-value - calc Z score

%addpath('/usr/local/matlab2008a/toolbox/matlab/specfun')
[h,pval,ci,zval] = ztest(max_ES_1,mean(EM_perms),std(EM_perms))

return
end
