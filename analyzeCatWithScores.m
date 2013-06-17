function analyzeCatWithScores(cat_ids, aspects, cat_pvalue_scores, regionLabels, htmlFileName,title, numOfGenesInCategory)


    

    GO = geneont('file', 'gene_ontology.obo.txt');
    cellularComponentCategories = strcmp('C', aspects); % not cellular components
    cat_ids = cat_ids(~cellularComponentCategories);

    sortedValues = nan(length(cat_ids), length(regionLabels));
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
       %fdrCorrectedValues = mafdr(relevantCategoriesScores, 'BHFDR', true);
       
       [sortedValues(:,i), sortIndices] = sort(relevantCategoriesScores);
%        [sortedValues(:,i), sortIndices] = sort(fdrCorrectedValues);
       sortedCatIndices(:,i) = cat_ids(sortIndices);
       sortedCatNames{i} = catNames(sortIndices);
       numOfgenesInCatAndRegion(:,i) = numOfGenesInCategory(sortIndices);
    end

    printCatScoreIntoHtml(regionLabels, sortedValues, sortedCatIndices, sortedCatNames, numOfgenesInCatAndRegion, title, htmlFileName)

end

