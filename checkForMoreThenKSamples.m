function [aboveK ,numOfSample] = checkForMoreThenKSamples(locationMatrix,experimentsSubjectMatrixLogical, ontologyObject,k)
    numOfPeople = size(experimentsSubjectMatrixLogical,2);
    [allChilds,~] = ontologyObject.allChildNodes();
    numOfSample = zeros(size(allChilds,1),numOfPeople);
    
    indicesToMean = allChilds';
    for i = 1:numOfPeople
        releventExperiments = experimentsSubjectMatrixLogical(:,i);
        experimentsInParentAreas = double( locationMatrix(releventExperiments,:) ) * double(indicesToMean);
        numOfSample(:,i) = sum(experimentsInParentAreas);
    end
%     for i=1:size(allChilds,1)
%         indicesToMean = allChilds(i,:)';
%         for numOfPersons = 1:length(locationMatrix)
%             releventExperiments = double(locationMatrix{numOfPersons}) * double(indicesToMean);
%             numOfSample(i,numOfPersons) = numOfSample(i) + sum(releventExperiments);
%         end
%     end
%     
    aboveK = all(numOfSample > k,2);
end