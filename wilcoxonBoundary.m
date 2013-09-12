function [rightBoundariesValues, leftBoundariesValues] = wilcoxonBoundary(X, rightBoundaries)
% X - rows are repetitions, columns are different variables
% returns the values at the position of the right and left boundaries for each variable
% return matrix has size #boundaries X #variables

    assert( all(rightBoundaries >= 0.5 &  rightBoundaries <= 1 ) );
    
    if any( any(isnan(X) ) )
        [rightBoundariesValues, leftBoundariesValues] = wilcoxonBoundaryWithNaN(X, rightBoundaries);
    else
        [rightBoundariesValues, leftBoundariesValues] = wilcoxonBoundaryWithoutNaN(X, rightBoundaries);
    end
    
end

function [rightBoundariesValues, leftBoundariesValues] = wilcoxonBoundaryWithoutNaN(X, rightBoundaries)
% X - rows are repetitions, columns are different variables
% returns the values at the position of the right and left boundaries for each variable
% return matrix has size #boundaries X #variables

    numberOfRepetitions = size(X,1);
    
    sortedX = sort(X, 1);
    
    leftBoundaries = 1 - rightBoundaries;
    
    rightBoundryIndecies = round(numberOfRepetitions * rightBoundaries); 
    leftBoundryIndecies = round(numberOfRepetitions * leftBoundaries); 
    
    rightBoundariesValues = sortedX(rightBoundryIndecies , :);
    leftBoundariesValues = sortedX(leftBoundryIndecies , :);
end

function [rightBoundariesValues, leftBoundariesValues] = wilcoxonBoundaryWithNaN(X, rightBoundaries)
% X - rows are repetitions, columns are different variables
% returns the values at the position of the right and left boundaries for each variable
% return matrix has size #boundaries X #variables

    numberOfVariables = size(X,2);
    numberOfBoundaries = length(rightBoundaries);

    rightBoundariesValues = nan(numberOfBoundaries, numberOfVariables);
    leftBoundariesValues = nan(numberOfBoundaries, numberOfVariables);
    
    leftBoundaries = 1 - rightBoundaries;
    
    for i =1 :numberOfVariables
        xWithoutNan = X(:,i);
        xWithoutNan = xWithoutNan( ~isnan(xWithoutNan));
        
        sorted_xi = sort(xWithoutNan);
        numberOfNotNanValues = length(xWithoutNan);

        rightBoundryIndecies = round(numberOfNotNanValues * rightBoundaries); 
        leftBoundryIndecies = round(numberOfNotNanValues * leftBoundaries); 

        rightBoundariesValues = sorted_xi(rightBoundryIndecies , :);
        leftBoundariesValues = sorted_xi(leftBoundryIndecies , :);
    end
end