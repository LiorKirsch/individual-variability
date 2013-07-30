function mainDistances()

%     load('easyFormatHumanData.mat');
% 
%     %distanceMatrixBetweenSamples = pdist2(experimentsDataMatrix,experimentsDataMatrix,'euclidean');
% 
%     expressionDistances = pdistAlternative(experimentsDataMatrix);
%     mniDistances = pdistAlternative(mni_xyz);
%     mriDistances = pdistAlternative(mri_voxel_xyz);
%       
%     save('dataWithDistance.mat');

    load('dataWithDistance.mat');
    load('humanOntologyObject.mat');
    load('subset.mat');

    numOfSamples = sum(experimentsLocationMatrix  ,1)';
   
    regionColors = humanOntology.structureColors(subsetNodesIndex,:);
    regionNames = humanOntology.structureLabels(subsetNodesIndex,4);



    numberOfPeople = size(experimentsSubjectMatrixLogical,2);
    numberOfRegions = size(experimentsLocationMatrix,2);
    
    [shortOntlogyDistances, longOntlogyDistances, meanDistances] = clacOntologyDistances(experimentsLocationMatrix, humanOntology);
    scatterDistances(expressionDistances, meanDistances, experimentsSubjectMatrixLogical, 'mean of the two distances to ontology common parent');
    scatterDistances(expressionDistances,  longOntlogyDistances, experimentsSubjectMatrixLogical, 'long distance to ontology common parent');
    scatterDistances(expressionDistances, shortOntlogyDistances, experimentsSubjectMatrixLogical, 'short distance to ontology common parent');
    
    
    scatterDistances(expressionDistances, mriDistances, experimentsSubjectMatrixLogical,'MRI based distance');
    scatterDistances(expressionDistances, mniDistances , experimentsSubjectMatrixLogical,'MNI based distance');
    
    
    
    regionandPeoplePrecision = nan(numberOfRegions,1);
    regionandPeoplePrecisionSTD = nan(numberOfRegions,1);

    for j = 1:numberOfRegions
        currentRegionExp =  experimentRegion(:,j);
        currentDistances = expressionDistances(currentRegionExp,:);
        currentDistances = currentDistances(:,currentRegionExp);
        currentExperimentsSubjectMatrixLogical = experimentsSubjectMatrixLogical(currentRegionExp,:);
        peoplePrecision = calcPrecisionAt1( currentDistances, currentExperimentsSubjectMatrixLogical );
        
        regionandPeoplePrecision(j) = mean(peoplePrecision);
        regionandPeoplePrecisionSTD(j) = std(peoplePrecision,0,2);
    end
        
    figure;
    hold on;
    bar(regionandPeoplePrecision);
    errorbar(regionandPeoplePrecision,regionandPeoplePrecisionSTD,'xr')
    
    xticklabel_rotate(1:length(regionandPeoplePrecision),40,regionNames, 'fontSize',14) ;
    
    figureHandle = gcf;
    set(findall(figureHandle,'type','text'),'fontSize',14);
    
    percentInCategory = calcPrecisionAt1( expressionDistances, experimentRegion );
    figure;
    bar(percentInCategory);
    xticklabel_rotate(1:length(percentInCategory),40,regionNames) ;

    plotDistances(experimentRegion, expressionDistances);
end

function [shortDistances, longDistances, meanDistances] = clacOntologyDistances(experimentsLocationMatrix, ontologyObj)
    numberOfRegions = size(experimentsLocationMatrix,2);
    numberOfSamples = size(experimentsLocationMatrix,1);
    
    shortDistances = nan(numberOfSamples);
    longDistances = nan(numberOfSamples);
    
    experimentRegionIndex = experimentsLocationMatrix * (1:size(experimentsLocationMatrix,2))';
    for i =1:numberOfSamples
        iIndex = experimentRegionIndex(i);
        shortDistances(i,:) = ontologyObj.shortDistanceToParent(iIndex, experimentRegionIndex);
        longDistances(i,:) = ontologyObj.longDistanceToParent(iIndex, experimentRegionIndex);
    end
    
    meanDistances = (shortDistances + longDistances) /2 ;
end
function plotDistances(experimentRegion, distanceMatrixBetweenSamples)
    experimentRegionIndex = experimentRegion * (1:size(experimentRegion,2))';
    [sortedValues,sortByRegion] = sort(experimentRegionIndex);
    sortedDistances = distanceMatrixBetweenSamples(sortByRegion, :);
    sortedDistances = sortedDistances(:,sortByRegion);
    figure;
    subplot(2,1,1);
    imagesc(sortedDistances);
    subplot(2,1,2);
    plot(sortedValues,'.');
    xlim([0, size(sortedDistances,1)]);
end

function distances = pdistAlternative(x)
    Qx=repmat(dot(x,x,2),1,size(x,1));
    distances = Qx+Qx'-2*(x*x');

    distances( logical(eye(size(distances)) )) = 0;
    distances = sqrt(distances);
end

function scatterDistances(expressionDistances, distances, experimentsSubjectMatrixLogical, label)
    numberOfSubjects = size(experimentsSubjectMatrixLogical,2);
    randI = randi(10000);
    
%     for i = 1:numberOfSubjects
%         releventExp = experimentsSubjectMatrixLogical(:,i);
%         releventExpressionDistances = expressionDistances(releventExp,releventExp);
%         releventDistances = distances(releventExp,releventExp);
%         pearsonCorr = corr(releventExpressionDistances(:),releventDistances(:),'type','Pearson');
%         spearmanCorr = corr(releventExpressionDistances(:),releventDistances(:),'type','Spearman');
%         
%         figure(randI + i);
%         scatter(releventExpressionDistances(:), releventDistances(:),'.');
%         xlabel('gene expression based distance');
%         ylabel(label);
%         title(sprintf('correlation: %g pearson, %g spearman', pearsonCorr, spearmanCorr));
%         
%     end
    
    spearmanCorr = corr(expressionDistances(:),distances(:),'type','Spearman');
    fprintf('%s \t %g \n', label, spearmanCorr);
end
    