function categoriesScores = clusteringScore(distanceBetweenCategories)
    % first sum over all genes
    distanctUsingAllFeatures = sum(distanceBetweenCategories,3);
    distanctUsingAllFeatures = squeeze( distanctUsingAllFeatures );
    
    distanceToItself  = diag(distanctUsingAllFeatures ) ;
    distanceToAllOthers  = sum(distanctUsingAllFeatures,2) - distanceToItself;
    categoriesScores = distanceToAllOthers./ distanceToItself;

end