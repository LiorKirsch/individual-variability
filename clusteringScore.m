function categoriesScores = clusteringScore(distanceBetweenCategories)
    % first sum over all genes
    distanctUsingAllFeatures = squeeze( sum(distanceBetweenCategories,3) );
    
    distanceToItself  = diag(distanctUsingAllFeatures ) ;
    distanceToAllOthers  = sum(distanctUsingAllFeatures,2) - distanceToItself;
    categoriesScores = distanceToAllOthers./ distanceToItself;

end