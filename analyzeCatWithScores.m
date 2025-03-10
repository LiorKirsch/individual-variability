function analyzeCatWithScores(cat_ids, aspects, cat_pvalue_scores, cat_score_seperate, regionLabels, go_gene_mat, geneNames,htmlFileName,title, numOfGenesInCategory, normalizedMeanGeneScore)

    GO = geneont('file', 'gene_ontology.obo.txt');
    cellularComponentCategories = strcmp('C', aspects); % not cellular components
    cat_ids = cat_ids(~cellularComponentCategories);

    sortedValues = nan(length(cat_ids), length(regionLabels));
    sortedValuesSeperate = nan(length(cat_ids), 6,length(regionLabels));
    sortedCatIndices = nan(length(cat_ids), length(regionLabels));
    numOfgenesInCatAndRegion = nan(length(cat_ids), length(regionLabels));
    sortedCatNames = cell(length(regionLabels),1);
    
    if ~exist('numOfGenesInCategory', 'var')
        numOfGenesInCategory = nan(size(cat_ids));
    end
    
    catNames = go_id2name(cat_ids, GO);
    
    for i =1:length(regionLabels)
       
       relevantCategoriesScores = cat_pvalue_scores{i};
       relevantCategoriesScores = relevantCategoriesScores(~cellularComponentCategories);
       
       relevantCategoriesScoresSeperate = cat_score_seperate{i};
       relevantCategoriesScoresSeperate = relevantCategoriesScoresSeperate(~cellularComponentCategories,:);
       
       %fdrCorrectedValues = mafdr(relevantCategoriesScores, 'BHFDR', true);
       
       [sortedValues(:,i), sortIndices] = sort(relevantCategoriesScores);
%        [sortedValues(:,i), sortIndices] = sort(fdrCorrectedValues);
        sortedValuesSeperate(:,:,i) = relevantCategoriesScoresSeperate(sortIndices,:);
       sortedCatIndices(:,i) = cat_ids(sortIndices);
       sortedCatNames{i} = catNames(sortIndices);
       sorted_gene_mat{i} = go_gene_mat(sortIndices,:);
       numOfgenesInCatAndRegion(:,i) = numOfGenesInCategory(sortIndices);
       
    end

    printCatScoreIntoHtml(regionLabels, sortedValues, sortedValuesSeperate, sortedCatIndices, sortedCatNames, sorted_gene_mat, geneNames, numOfgenesInCatAndRegion, title, htmlFileName, normalizedMeanGeneScore)

end

