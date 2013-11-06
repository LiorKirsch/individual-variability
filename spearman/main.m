addpath('/home/lab/gal/develop/matlab');
%
% Analayze global gradients in human and mouse 
%

init 

datanames = {'kang','mouse', 'human6'}
num_ds = length(datanames);

for i_ds = 1:num_ds
    [spearman_scores, scores_rnd, gene_symbols] = compute_spearman(datanames{i_ds}, parms);

    figure(i_ds); clf; hold on;
    bins = [-1:0.01:1];
    a = histc(spearman_scores, bins);
    a = a / sum(a);
    plot(bins, a, 'Color', [0.5 0.5 0.99], 'LineWidth', 3);
    ax = axis;
    
    a = histc(cell2mat(scores_rnd'), bins);
    a = a / sum(a);    
    plot(bins , a, 'Color', [0.7 0.2 0.2], 'LineWidth', 3);
    axis(ax)
    legend(datanames{i_ds} , 'random');
end






