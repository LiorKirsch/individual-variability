function [aboveK ,numOfSample] = checkForMoreThenKSamples(locationMatrix, ontologyObject,k)
    [allChilds,~] = ontologyObject.allChildNodes();
    numOfSample = zeros(size(allChilds,1),1);
    
        
    for i=1:size(allChilds,1)
        indicesToMean = allChilds(i,:)';
        for numOfPersons = length(locationMatrix)
            releventExperiments = double(locationMatrix{numOfPersons}) * double(indicesToMean);
            numOfSample(i) = numOfSample(i) + sum(releventExperiments);
        end
    end
    
    aboveK = numOfSample >= k;
end