function [linearRegressionWeights, lassoWeights, ordinalWeights] = findWeightsOfFit(X, y)

    numberOfSamples = size(X,1);
    numberOfFeatures = size(X,2);
    
    assert( length(y) == numberOfSamples);
    
    linearRegressionWeights = polyfit(X, y, 1);
    linearRegressionWeights2 = [ones(numberOfSamples, 1), X] \ y;
    
    
    [lassoWeights, FitInfo] = lasso(X,y,'CV',10);
    lassoPlot(B,FitInfo,'PlotType','CV');
    
    ordinalWeights = mnrfit(X,y, 'model','ordinal');
    
end