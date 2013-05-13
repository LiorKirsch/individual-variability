function [anovaMatrix, samplesInRegion]= calcAnovaInStructures(experimentData, locationMatrix, ontologyObject, subsetNodesIndex)

    numOfSubjects = length(experimentData);
    assert(length(experimentData) == length(locationMatrix));
    for i = 1:length(experimentData)
        assert(size(experimentData{i}.expressionMatrix,2) == size(locationMatrix{i},1));
        assert(size(locationMatrix{i},2) == size(ontologyObject.dependencyMatrix,1) );
    end
    
    
    [allChilds,~] = ontologyObject.allChildNodes();
    childsOfSelectedNodes = allChilds(subsetNodesIndex,:);
    anovaMatrix = nan(size(experimentData{1}.expressionMatrix,1), size(childsOfSelectedNodes,1) );
    samplesInRegion = nan( size(childsOfSelectedNodes,1), numOfSubjects);
    
    for i=1:size(childsOfSelectedNodes,1)
%         indicesToMean = false(size(locationMatrix,2),1);
%         indicesToMean(i) = true;
        indicesToMean = childsOfSelectedNodes(i,:)';
        [anovaMatrix(:,i), samplesInRegion(i,:)] = meanExpressionOfIndices(experimentData, locationMatrix, indicesToMean);
    end
end

function [anovaScoresForGenes, samplesInRegion]= meanExpressionOfIndices(experimentData, locationMatrix, indicesToMean)
    
    for i = 1:length(locationMatrix)
        assert(size(locationMatrix{i},2) == size(indicesToMean,1));
    end
    
    numOfsubjects = length(experimentData);
    numOfGenes = size(experimentData{1}.expressionMatrix,1);
    releventExprimentsData  =nan(numOfGenes,0);
    experimentToSubjectID = [];
    samplesInRegion = nan(length(experimentData),1);
    
    for i = 1:numOfsubjects
        releventExperiments = double(locationMatrix{i}) * double(indicesToMean);
        releventExperiments = logical(releventExperiments);
        currentSubjectExperiments = experimentData{i}.expressionMatrix(:,releventExperiments);
        
        samplesInRegion(i) = size(currentSubjectExperiments,2);
        releventExprimentsData = cat(2,releventExprimentsData,currentSubjectExperiments);
        currentGroup =  ones( size(currentSubjectExperiments, 2),1) * i;
        experimentToSubjectID = cat(1,experimentToSubjectID,currentGroup);
    end
    
    
    anovaScoresForGenes = computeAnova(releventExprimentsData, experimentToSubjectID);
end

function anovaScoresForGenes = computeAnova(releventExprimentsData, groupData)
    numOfGenes = size(releventExprimentsData,1);
    anovaScoresForGenes = nan(numOfGenes,1);
    parfor i =1:numOfGenes %number of genes
        anovaScoresForGenes(i) = anova1(releventExprimentsData(i,:), groupData,'off') ;
    end
end