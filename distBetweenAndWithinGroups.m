function distanceBetweenCategories = distBetweenAndWithinGroups(dataMatrix, samplesInGroup)
% for each feature calculate the distance between categories (subjects)    
    numberOfCategories = size(samplesInGroup,2);
    numberOfFeatures = size(dataMatrix,2);
    
    distanceBetweenCategories = zeros(numberOfCategories, numberOfCategories, numberOfFeatures);
    
    for i =1:numberOfCategories
        for j =i:numberOfCategories
            distanceBetweenCategories(i, j ,:) = calcDistBetweenGroups(dataMatrix, samplesInGroup,i, j);
            distanceBetweenCategories(j, i ,:) = calcDistBetweenGroups(dataMatrix, samplesInGroup,i, j);
        end
    end
end


function score = calcDistBetweenGroups(dataMatrix, samplesCategory, firstCatIndex, secondCatIndex)

    group1Indcies = samplesCategory(:,firstCatIndex);
    group2Indcies = samplesCategory(:,secondCatIndex);
    
    group1Size = sum(group1Indcies);
    group2Size = sum(group2Indcies);
    
    group1data = dataMatrix(group1Indcies,:);
    group2data = dataMatrix(group2Indcies,:);
    group1Var = var(group1data,1,1);
    group2Var = var(group2data,1,1);
    group1Mean = mean(group1data,1);
    group2Mean = mean(group2data,1);
    
    score = group1Var + group2Var + (group1Mean - group2Mean).^2;
    
end