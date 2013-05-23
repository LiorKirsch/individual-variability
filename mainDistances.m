function mainDistances()
%     load('pcaData.mat', 'experimentRegion', 'experimentsDataMatrix','experimentsLocationMatrix','experimentsSubjectMatrix','regionNames');
%     distanceMatrixBetweenSamples = pdist2(experimentsDataMatrix,experimentsDataMatrix,'euclidean');
%     save('dataWithDistance.mat');

    load('dataWithDistance.mat');
    
    numberOfPeople = 6;
    numberOfRegions = size(experimentRegion,2);
    
    experimentsSubjectMatrixLogical = false(size(experimentsSubjectMatrix,1), numberOfPeople);
    for i =1:numberOfPeople
        experimentsSubjectMatrixLogical(:,i) = experimentsSubjectMatrix == i;
    end
    
    
    regionandPeoplePrecision = nan(numberOfRegions,1);
    regionandPeoplePrecisionSTD = nan(numberOfRegions,1);

    for j = 1:numberOfRegions
        currentRegionExp =  experimentRegion(:,j);
        currentDistances = distanceMatrixBetweenSamples(currentRegionExp,:);
        currentDistances = currentDistances(:,currentRegionExp);
        currentExperimentsSubjectMatrixLogical = experimentsSubjectMatrixLogical(currentRegionExp,:);
        peoplePrecision = showPrecisionAt1( currentDistances, currentExperimentsSubjectMatrixLogical ,1:6);
        
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
    
    percentInCategory = showPrecisionAt1( distanceMatrixBetweenSamples, experimentRegion ,1:15);
    figure;
    bar(percentInCategory);
    xticklabel_rotate(1:length(percentInCategory),40,regionNames) ;

    plotDistances(experimentRegion, distanceMatrixBetweenSamples);
end

function percentInCategory = showPrecisionAt1( distanceMatrixBetweenSamples, experimentCategories ,groupLabels)

    distanceMatrixBetweenSamples(1:size(distanceMatrixBetweenSamples,1)+1:end) = inf;
    [minValues, minIndices] = min(distanceMatrixBetweenSamples,[],1);

    nearestNeighborGroup = experimentCategories(minIndices,:);
    numberOfCategories = size(experimentCategories,2);
    
    samplesWhoseNearestNeighborIsInTheSameCateogry = false(size(experimentCategories));
    for i =1:numberOfCategories
        samplesInRegion = experimentCategories(:,i);
        nearestNeighborInRegion = nearestNeighborGroup(:,i);
        samplesWhoseNearestNeighborIsInTheSameCateogry(:,i) = samplesInRegion & nearestNeighborInRegion;
    end

    percentInCategory = sum(samplesWhoseNearestNeighborIsInTheSameCateogry,1) ./ sum(experimentCategories,1);
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