function mainDistances()
%     load('pcaData.mat', 'experimentRegion', 'experimentsDataMatrix','experimentsLocationMatrix','experimentsSubjectMatrix','regionNames','mni_xyz', 'mri_voxel_xyz');
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
    aboveK = numOfSamples >= 10;
    
    subsetNodesIndex = subsetNodesIndex & aboveK;
    %subsetNodesIndex = true(size(aboveK));
    regionColors = humanOntology.structureColors(subsetNodesIndex,:);
    regionNames = humanOntology.structureLabels(subsetNodesIndex,4);



    numberOfPeople = 6;
    numberOfRegions = size(experimentRegion,2);
    
    experimentsSubjectMatrixLogical = false(size(experimentsSubjectMatrix,1), numberOfPeople);
    for i =1:numberOfPeople
        experimentsSubjectMatrixLogical(:,i) = experimentsSubjectMatrix == i;
    end
    
    scatterDistances(expressionDistances, mniDistances , mriDistances, experimentsSubjectMatrixLogical);
    
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
    distances = sqrt(Qx+Qx'-2*(x*x'));
end

function scatterDistances(expressionDistances, mniDistances , mriDistances, experimentsSubjectMatrixLogical)
    numberOfSubjects = size(experimentsSubjectMatrixLogical,2);
    for i = 1:numberOfSubjects
        releventExp = experimentsSubjectMatrixLogical(:,i);
        releventExpressionDistances = expressionDistances(releventExp,releventExp);
        releventMRIDistances = mriDistances(releventExp,releventExp);
        releventMNIDistances = mniDistances(releventExp,releventExp);
        
        figure(i*10);
        scatter(releventExpressionDistances(:), releventMNIDistances(:),'.');
        xlabel('gene expression based distance');
        ylabel('MNI based distance');

        figure(i*20);
        scatter(releventExpressionDistances(:), releventMRIDistances(:),'.');
        xlabel('gene expression based distance');
        ylabel('MRI based distance');
    end
end
    