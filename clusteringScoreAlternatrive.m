function [categoriesScores, categoriesScoresUsingFeatures] = clusteringScoreAlternatrive(distanceBetweenCategories,fullDistances)
    % first sum over all genes
    distanctUsingSelectedFeatures = sum(distanceBetweenCategories,3);
    distanctUsingSelectedFeatures = squeeze( distanctUsingSelectedFeatures );
    distanctUsingAllFeatures = fullDistances - distanctUsingSelectedFeatures;
    
    categoriesScores = calcCategoriesScores(distanctUsingAllFeatures);
    categoriesScoresUsingFeatures = calcCategoriesScores(distanctUsingSelectedFeatures);
end

function categoriesScores = calcCategoriesScores(distancesUsingFeatures)
    distanceToItself  = diag(distancesUsingFeatures ) ;
    distanceToAllOthers  = sum(distancesUsingFeatures,2) - distanceToItself;
    categoriesScores = distanceToAllOthers./ distanceToItself;
end