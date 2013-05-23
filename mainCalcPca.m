function mainCalcPca()

% 
% load('expressionRawDataWithOntology.mat', 'expressionDataSelected',  'humanOntology','selectedProbesData','locationMatrix');
% load('humanOntologyObject.mat');
% load('subset.mat');
% 
% [aboveK ,numOfSamples] = checkForMoreThenKSamples(locationMatrix, humanOntology,10);
% subsetNodesIndex = subsetNodesIndex & aboveK;
% %subsetNodesIndex = true(size(aboveK));
% regionColors = humanOntology.structureColors(subsetNodesIndex,:);
% regionNames = humanOntology.structureLabels(subsetNodesIndex,4);
% 
% [experimentsDataMatrix, experimentsLocationMatrix, experimentsSubjectMatrix, probeCode, experimentRegion] = calcPCAinStructures(expressionDataSelected, locationMatrix, humanOntology,subsetNodesIndex);
% [coeff,score,latent] = pca(experimentsDataMatrix);
% save('pcaData.mat','experimentsDataMatrix', 'experimentsLocationMatrix', 'experimentsSubjectMatrix', 'probeCode', 'experimentRegion', 'coeff','score', 'latent','regionNames','regionColors');
% %save('pcaDataAllRegions.mat','experimentsDataMatrix', 'experimentsLocationMatrix', 'experimentsSubjectMatrix', 'probeCode', 'experimentRegion', 'coeff','score', 'latent','regionNames','regionColors');

load('pcaData.mat');
%load('pcaDataAllRegions.mat');
numOfSubjects = max(max(experimentsSubjectMatrix));

regionColors = regionColors/255;
regionColors = createColorMap(length(regionNames));

drawPCAregionColor(score, experimentRegion, experimentsSubjectMatrix, numOfSubjects, regionColors,regionNames);
drawPCA(score, experimentRegion, experimentsSubjectMatrix, numOfSubjects,regionColors, regionNames);

drawPCApeopleColor(score, experimentsSubjectMatrix,0,'',latent);
clacPCAforEachAreaSeperately(experimentsDataMatrix, experimentRegion, experimentsSubjectMatrix, regionNames);
colorbar
end

function clacPCAforEachAreaSeperately(experimentsDataMatrix, experimentRegion, experimentsSubjectMatrix, regionNames)
    numberOfAreas = size(experimentRegion,2);
    for i = 1:numberOfAreas
       releventExperementsIncdices =  experimentRegion(:,i);
       releventExperements = experimentsDataMatrix(releventExperementsIncdices,:);
       releventExperimentsSubjectMatrix = experimentsSubjectMatrix(releventExperementsIncdices);
       [coeff,score,latent] = pca(releventExperements);
       
       drawPCApeopleColor(score, releventExperimentsSubjectMatrix,i*10,regionNames{i},latent);
    end
end

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

function [experimentsDataMatrix, experimentsLocationMatrix, experimentsSubjectMatrix, probeCode, experimentRegion]= calcPCAinStructures(experimentData, locationMatrix, ontologyObject, subsetNodesIndex)

    numOfSubjects = length(experimentData);
    assert(length(experimentData) == length(locationMatrix));
    for i = 1:length(experimentData)
        assert(size(experimentData{i}.expressionMatrix,2) == size(locationMatrix{i},1));
        assert(size(locationMatrix{i},2) == size(ontologyObject.dependencyMatrix,1) );
    end
    
    
    [allChilds,~] = ontologyObject.allChildNodes();
    childsOfSelectedNodes = allChilds(subsetNodesIndex,:);
    
    [experimentsDataMatrix, experimentsLocationMatrix, experimentsSubjectMatrix, probeCode] = joinDataToOneBigMatrix(experimentData, locationMatrix);
    experimentsDataMatrix = experimentsDataMatrix';
    
    experimentRegion =  false(size(experimentsDataMatrix,1), size(childsOfSelectedNodes,1));
    for i=1:size(childsOfSelectedNodes,1) % loop over the high category areas
        indicesToMean = childsOfSelectedNodes(i,:)';
        releventExperiments = double(experimentsLocationMatrix) * double(indicesToMean);
        experimentRegion(:,i) = logical(releventExperiments);
    end
    
    sampleWithNoLocation = ~any(experimentRegion,2);
     experimentsDataMatrix = experimentsDataMatrix(~sampleWithNoLocation,:);
     experimentsLocationMatrix = experimentsLocationMatrix(~sampleWithNoLocation,:);
     experimentsSubjectMatrix = experimentsSubjectMatrix(~sampleWithNoLocation,:);
     experimentRegion = experimentRegion(~sampleWithNoLocation,:);
    
end

function drawPCAregionColor(score, experimentRegion, experimentsSubjectMatrix, numOfSubjects, regionColors,regionNames)
    numberOfRegions = size(experimentRegion,2);
    experimentRegionIndex = experimentRegion *  (1:numberOfRegions)';
    experimentRegionLabel = regionNames(experimentRegionIndex,:);
    experimentRegionColor = regionColors(experimentRegionIndex,:);
    markerType = {'o' ,'^','v','x','s','d', '>','+','p' ,'h','s' ,'<','x','o','d'};

    
    for j=1: numberOfRegions
       tmpIndex = experimentRegion(:,j);
       tmpRegion =  experimentRegionLabel(tmpIndex);
       tmpColor = experimentRegionColor(tmpIndex,:);
       tmpIExperiments = score(tmpIndex,:);
       figure(1);
       hold on;
       assert( all(tmpColor(1,:) == regionColors(j,:) ) );
       scatter(tmpIExperiments(:,1),tmpIExperiments(:,2),10,tmpColor(1,:),markerType{j});

        figure(2);
        hold on;
        scatter(tmpIExperiments(:,2),tmpIExperiments(:,3),10,tmpColor(1,:),markerType{j});
        figure(3);
        hold on;
        scatter(tmpIExperiments(:,3),tmpIExperiments(:,4),10,tmpColor(1,:),markerType{j});

    end
    figure(1);
    xlabel('Principle component #1');     ylabel('Principle component #2');
    legend(regionNames);
    figureHandle = gcf;     set(findall(figureHandle,'type','text'),'fontSize',14);
    figure(2);
    xlabel('Principle component #2');     ylabel('Principle component #3');
    legend(regionNames);
    figureHandle = gcf;     set(findall(figureHandle,'type','text'),'fontSize',14);
    figure(3);
    xlabel('Principle component #3');     ylabel('Principle component #4');
    legend(regionNames);
    figureHandle = gcf;     set(findall(figureHandle,'type','text'),'fontSize',14);
end

function colors = createColorMap(numOfElements)
    colormapValues = colormap(jet);
    indexOfColor = round(1:size(colormapValues,1)/numOfElements:size(colormapValues,1));
    colors = colormapValues(indexOfColor,:);
end
function drawPCA(score, experimentRegion, experimentsSubjectMatrix, numOfSubjects, regionColors,regionNames)
    numberOfRegions = size(experimentRegion,2);
    markerType = {'+', 'o', '*' , 'x' , 's', 'd'};
    experimentRegionIndex = experimentRegion *  (1:numberOfRegions)';
    experimentRegionLabel = regionNames(experimentRegionIndex,:);
    experimentRegionColor = regionColors(experimentRegionIndex,:);
    experimentRegionColor = experimentRegionColor/255;
    for i =1:numOfSubjects
        subjectIExperimentsIndices = experimentsSubjectMatrix == i ;
        
        subjectIExperiments = score(subjectIExperimentsIndices,:);
        subjectIAreas = experimentRegionIndex(subjectIExperimentsIndices) ;
        subjectIRegionNames = experimentRegionLabel(subjectIExperimentsIndices) ;
        subjectIColors = experimentRegionColor(subjectIExperimentsIndices,:) ;
         
        figure(1);
        hold on;
        scatter(subjectIExperiments(:,1),subjectIExperiments(:,2),60,subjectIAreas,markerType{i});
        xlabel('Principle component #1');    ylabel('Principle component #2');
        figure(2);
        hold on;
        scatter(subjectIExperiments(:,2),subjectIExperiments(:,3),60,subjectIAreas,markerType{i});
        xlabel('Principle component #2');    ylabel('Principle component #3');
        figure(3);
        hold on;
        scatter(subjectIExperiments(:,3),subjectIExperiments(:,4),60,subjectIAreas,markerType{i});
        xlabel('Principle component #3');    ylabel('Principle component #4');
    end

end

function drawPCApeopleColor(score, experimentsSubjectMatrix,figureIndex,titleString, latent)
    
    if ~exist('figureIndex','var')
        figureIndex = 0;
    end
    figure(figureIndex+ 1);
    subplot(2,2,1);
    scatter(score(:,1),score(:,2),50,experimentsSubjectMatrix,'filled');
    xlabel('Principle component #1');    ylabel('Principle component #2');
    subplot(2,2,2);
    scatter(score(:,2),score(:,3),50,experimentsSubjectMatrix,'filled');
    xlabel('Principle component #2');    ylabel('Principle component #3');
    subplot(2,2,3);
    scatter(score(:,3),score(:,4),50,experimentsSubjectMatrix,'filled');
    xlabel('Principle component #3');    ylabel('Principle component #4');
    
    subplot(2,2,4);
    plot(cumsum(latent)/sum(latent),'.-')
    title('variance explained');
    xlabel('eigen value index');
    xlim([-1 60]);
%     figure(figureIndex + 4);
%     scatter(score(:,1),score(:,3),50,experimentsSubjectMatrix,'filled');
%     title(titleString, ' - pca 12']);
             
   
end

function [experimentsDataMatrix, experimentsLocationMatrix, experimentsSubjectMatrix, probeCode] = joinDataToOneBigMatrix(experimentData, locationMatrix)
    numberOfProbes = size(experimentData{1}.expressionMatrix, 1);
    numberOfSamples = nan(length(experimentData),1);
    numberOfSubjects = length(experimentData);
    numberOfAreas = size(locationMatrix{1}, 2);
    
    for i=1:numberOfSubjects
        numberOfSamples(i) = size(experimentData{i}.expressionMatrix, 2);
    end
    totalNumberOfSamples = sum(numberOfSamples);
    
    experimentsDataMatrix = nan(numberOfProbes, totalNumberOfSamples);
    experimentsLocationMatrix = nan(totalNumberOfSamples, numberOfAreas);
    experimentsSubjectMatrix = nan(totalNumberOfSamples,1);

    lastIndex = 1;
    probeCode = experimentData{1}.probeCode;
    for i = 1:numberOfSubjects
       assert(all(probeCode == experimentData{i}.probeCode));
       indicesOfIntrests = lastIndex: lastIndex + numberOfSamples(i) - 1;
       experimentsDataMatrix(:, indicesOfIntrests) = experimentData{i}.expressionMatrix;
       experimentsLocationMatrix( indicesOfIntrests , :) = locationMatrix{i};
       experimentsSubjectMatrix ( indicesOfIntrests) =  i ;
       lastIndex = lastIndex + numberOfSamples(i);
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
    someViewingForAnova(anovaScoresForGenes, releventExprimentsData , groupData);
end

