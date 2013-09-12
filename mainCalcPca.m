function mainCalcPca()


% load('easyFormatHumanData.mat');
% load('humanOntologyObject.mat');
% load('subset.mat');
% 
% [aboveK ,numOfSamples] = checkForMoreThenKSamples(experimentsLocationMatrix, experimentsSubjectMatrixLogical, humanOntology,5);
% subsetNodesIndex = subsetNodesIndex & aboveK;
% %subsetNodesIndex = true(size(aboveK));
% regionColors = humanOntology.structureColors(subsetNodesIndex,:);
% regionNames = humanOntology.structureLabels(subsetNodesIndex,4);
% 
% [experimentsDataMatrix, experimentsLocationMatrix, experimentsSubjectMatrix, experimentRegion, mni_xyz, mri_voxel_xyz] = groupByRegions(experimentsDataMatrix, experimentsLocationMatrix, experimentsSubjectMatrixLogical,humanOntology,subsetNodesIndex, mni_xyz, mri_voxel_xyz);
% [coeff,score,latent] = pca(experimentsDataMatrix);
% save('pcaData.mat','experimentsDataMatrix', 'experimentsLocationMatrix', 'experimentsSubjectMatrix', 'selectedProbesData', 'experimentRegion', 'coeff','score', 'latent','regionNames','regionColors','mni_xyz', 'mri_voxel_xyz');
% %save('pcaDataAllRegions.mat','experimentsDataMatrix', 'experimentsLocationMatrix', 'experimentsSubjectMatrix', 'selectedProbesData', 'experimentRegion', 'coeff','score', 'latent','regionNames','regionColors','mni_xyz', 'mri_voxel_xyz');
%%
%load('pcaData.mat');
load('pcaSparseData.mat');

load('humanOntologyObject.mat');
%load('pcaDataAllRegions.mat');
numOfSubjects = size(experimentsSubjectMatrix,2);

regionColors = regionColors/255;
regionColors = createColorMap(length(regionNames));

regionColors = humanOntology.getColorByRegionName(regionNames);

%plotMds(experimentsDataMatrix,experimentRegion, experimentsSubjectMatrix, regionNames,regionColors);

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

function [experimentsDataMatrix, locationMatrix, experimentsSubjectMatrixLogical, experimentRegion, mni_xyz, mri_voxel_xyz]= groupByRegions(experimentsDataMatrix, locationMatrix, experimentsSubjectMatrixLogical,ontologyObject, subsetNodesIndex, mni_xyz, mri_voxel_xyz)

    assert(size(experimentsDataMatrix,1) == size(locationMatrix,1));
    assert(size(locationMatrix,2) == size(ontologyObject.dependencyMatrix,1) );
    
    experimentRegion = getRegionOfGrossStructures(subsetNodesIndex, locationMatrix, ontologyObject);

    sampleWithNoLocation = ~any(experimentRegion,2);
    experimentsDataMatrix = experimentsDataMatrix(~sampleWithNoLocation,:);
    locationMatrix = locationMatrix(~sampleWithNoLocation,:);
    experimentsSubjectMatrixLogical = experimentsSubjectMatrixLogical(~sampleWithNoLocation,:);
    experimentRegion = experimentRegion(~sampleWithNoLocation,:);
    mni_xyz = mni_xyz(~sampleWithNoLocation,:);
    mri_voxel_xyz = mri_voxel_xyz(~sampleWithNoLocation,:);
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
    legend(regionNames,'Location','SouthWest');
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

function drawPCA(score, experimentRegion, experimentsSubjectMatrix, numOfSubjects, regionColors,regionNames)
    numberOfRegions = size(experimentRegion,2);
    markerType = {'+', 'o', '*' , 'x' , 's', 'd'};
    experimentRegionIndex = experimentRegion *  (1:numberOfRegions)';
    experimentRegionLabel = regionNames(experimentRegionIndex,:);
    experimentRegionColor = regionColors(experimentRegionIndex,:);
    experimentRegionColor = experimentRegionColor/255;
    for i =1:numOfSubjects
        subjectIExperimentsIndices = experimentsSubjectMatrix(:,i) ;
        
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
    experimentsSubjectMatrix = experimentsSubjectMatrix * (1:size(experimentsSubjectMatrix,2))';
    figure(figureIndex+ 1);
    subplot(2,2,1);
    scatter(score(:,1),score(:,2),50,experimentsSubjectMatrix);
    xlabel('Principle component #1');    ylabel('Principle component #2');
    subplot(2,2,2);
    scatter(score(:,2),score(:,3),50,experimentsSubjectMatrix);
    xlabel('Principle component #2');    ylabel('Principle component #3');
    subplot(2,2,3);
    scatter(score(:,3),score(:,4),50,experimentsSubjectMatrix);
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

