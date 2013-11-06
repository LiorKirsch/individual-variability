function [scores, scores_rnd, gene_symbols] = compute_spearman(dataname, params)
% 

    [expression, gross_region_vec, gene_symbols] = load_expression_and_regions(dataname, params);

    scores = corr(expression, gross_region_vec,'type','Spearman');
    
    [sortedScores ,sortInd] = sort(scores, 'descend');
    writeCellToFile([dataname,'_spearman.txt'], gene_symbols(sortInd) );
    
    
    [sortedScores ,sortInd] = sort(abs(scores), 'descend');
    writeCellToFile([dataname,'_abs_spearman.txt'], gene_symbols(sortInd) );
    

    % Randomize            
    old_seed = rand('state'); %#ok
    for i =1:3
        rand('state', i);
        gross_region_vec_rnd = randrows(gross_region_vec);
        scores_rnd{i} = corr(expression, gross_region_vec_rnd,'type','Spearman');
    end
    
    rand('state', old_seed); %#ok    
end


