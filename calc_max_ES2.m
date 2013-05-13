function [max_ES,ES] = calc_max_ES2(sorted_pred,sorted_labels,p)

% Lior Kirsch based on Ossnat mar. 2011
% sorted_pred - the rank of the prediction (r_i in the paper)
%               all values are possitive - no need in abs()
% sorted_labels - indicator vector {0,1} wheter the gene belongs to the
% category 
% p - should be 1 or 0, 1 normalizing by the sum, 0 standard kolmogorov
% smirnov statiastics 
% NH = # gene in Set S
% N = #genes

%
% Paper:"Gene set enrichment analysis: A knowledge-based approach for
% interpreting genome-wide expression profiles"
% http://www.pnas.org/content/102/43/15545.long

switch(p)
    case 0,
        S_i = find(sorted_labels>0);
        N_R = length(S_i);
        N_R_vec = N_R*ones(size(sorted_pred));
        N_minus_N_H = length(sorted_pred)-length(S_i); % total # of genes - genes in S
        r_j_p = sorted_labels; 
        P_hit = cumsum(r_j_p)./N_R_vec;

        P_miss_cont = (1/N_minus_N_H)*ones(size(sorted_pred)).*(~sorted_labels);
        P_miss =  cumsum(P_miss_cont);

    case 1,
        S_i = find(sorted_labels>0);
        N_R = sum(sorted_pred(S_i));
%         N_R_vec = N_R*ones(size(sorted_pred));
        N_minus_N_H = length(sorted_pred)-length(S_i); % total # of genes - genes in S
        r_j_p = sorted_labels.*sorted_pred;
%         r_j_p = zeros(size(sorted_pred));
%         r_j_p(sorted_labels) = sorted_pred(sorted_labels);
        P_hit = cumsum(r_j_p)/N_R;
        P_miss_cont = (1/N_minus_N_H)*(~sorted_labels);
        P_miss =  cumsum(P_miss_cont);
end

ES = P_hit-P_miss;
max_ES = max(ES);


return
end
