function mainAgeVaiability()
    load('easyFormatHumanData_withAges.mat');
    load('subset.mat');
    load('humanOntologyObject.mat');
    randomAgeScores = load('randomSizeAgePrecision.mat');
    regionSampleTreshold = 50;
    
    regionColors = humanOntology.structureColors(subsetNodesIndex,:);
    regionNames = humanOntology.structureLabels(subsetNodesIndex,4);
    experimentRegion = getRegionOfGrossStructures(subsetNodesIndex, experimentsLocationMatrix, humanOntology);

    experimentsInRegion = sum(experimentRegion,1);
    aboveTreshold = (experimentsInRegion > regionSampleTreshold);
    experimentRegion = experimentRegion(:,aboveTreshold);
    regionNames = regionNames(aboveTreshold);
    regionColors = regionColors(aboveTreshold,:);
    experimentsInRegion = experimentsInRegion(aboveTreshold);

    numberOfPeople = size(experimentsSubjectMatrixLogical,2);
    numberOfSamples = size(experimentsSubjectMatrixLogical,1);
    numberOfRegions = size(experimentRegion,2);

    assert(  numberOfSamples == size(experimentsDataMatrix,1) );

    allSamples = true(size(experimentRegion,1),1);
    
    experimentsDataMatrixNormalized = normalizeMatrixByDim(experimentsDataMatrix, 1);

    
    %[categoriesScoresAllSamples, distanceBetweenGroupsAllSamples] = calcRegionBetweenSubjectDist(experimentsSubjectMatrixLogical, experimentsDataMatrix, allSamples);
    [categoriesScoresAllSamples, distanceBetweenGroupsAllSamples] = calcRegionBetweenSubjectDist(experimentsSubjectMatrixLogical, experimentsDataMatrixNormalized, allSamples);
    numberOfEra = size(agesDescription,1);

    eraDistanceScoreAllSamples = nan(numberOfEra,numberOfPeople);
    onlyEraDistanceScoreAllSamples = nan(numberOfEra,numberOfPeople);
    
    randomScoresMeanPeople = squeeze(mean(randomAgeScores.percentInCategoryUsingGenes,2));
    randomScoresMeanPeopleWithoutGenes = squeeze(mean(randomAgeScores.percentInCategory,2));
    
    for i = 1:numberOfEra

        eraName = agesDescription{i}; 

        genesInEra = logical(ages(:,i));
        currentDistanceBetweenGroups = distanceBetweenGroupsAllSamples(:,:,:,~genesInEra);
        currentDistanceBetweenGroupsOnlyEra = distanceBetweenGroupsAllSamples(:,:,:,genesInEra);

        currentAreaDistanceBetweenGroups = squeeze(currentDistanceBetweenGroups(1,:,:,:));
        eraDistanceScoreAllSamples(i,:) = clusteringScore( currentAreaDistanceBetweenGroups );

        currentAreaOnlyCatDistanceBetweenGroups = squeeze(currentDistanceBetweenGroupsOnlyEra(1,:,:,:));
        onlyEraDistanceScoreAllSamples(i,:) = clusteringScore( currentAreaOnlyCatDistanceBetweenGroups );

        printPercentCounter(i, numberOfEra );
    end

    %[categoriesScores, distanceBetweenGroups] = calcRegionBetweenSubjectDist(experimentsSubjectMatrixLogical, experimentsDataMatrix, experimentRegion);
    [categoriesScores, distanceBetweenGroups] = calcRegionBetweenSubjectDist(experimentsSubjectMatrixLogical, experimentsDataMatrixNormalized, experimentRegion);

    
    eraDistanceScore = nan(numberOfEra,numberOfPeople,numberOfRegions);
    onlyEraDistanceScore = nan(numberOfEra,numberOfPeople,numberOfRegions);

    for i = 1:numberOfEra

        eraName = agesDescription{i}; 

        genesInEra = logical(ages(:,i));
        currentDistanceBetweenGroups = distanceBetweenGroups(:,:,:,~genesInEra);
        currentDistanceBetweenGroupsOnlyEra = distanceBetweenGroups(:,:,:,genesInEra);

        for m =1:numberOfRegions
            currentAreaDistanceBetweenGroups = squeeze(currentDistanceBetweenGroups(m,:,:,:));
            eraDistanceScore(i,:,m) = clusteringScore( currentAreaDistanceBetweenGroups );

            currentAreaOnlyCatDistanceBetweenGroups = squeeze(currentDistanceBetweenGroupsOnlyEra(m,:,:,:));
            onlyEraDistanceScore(i,:,m) = clusteringScore( currentAreaOnlyCatDistanceBetweenGroups );
        end
        printPercentCounter(i, numberOfEra );
    end

    drawVariability(eraDistanceScoreAllSamples , onlyEraDistanceScoreAllSamples, 'whole brain','gene era');

    for m = 1:numberOfRegions
       figure(m);
       areaRandomScores = squeeze(randomScoresMeanPeople(:,m ,:));
       areaRandomScores = areaRandomScores';
       areaRandomScoresWithout = squeeze(randomScoresMeanPeopleWithoutGenes(:,m ,:));
       areaRandomScoresWithout = areaRandomScoresWithout';
       
       currentTitle = sprintf('%s (%d genes)',regionNames{m}, experimentsInRegion(m));
       %drawVariability(eraDistanceScore(:,:,m) , onlyEraDistanceScore(:,:,m), '','gene era');
       drawVariability(areaRandomScoresWithout , areaRandomScores, '','gene era');
       hold on;
       plot( mean(onlyEraDistanceScore(:,:,m),2) );
       annotation('textbox', [0 0.9 1 0.1], 'String', currentTitle, 'EdgeColor', 'none', 'HorizontalAlignment', 'center','fontsize',20);
    end
    
    
%     for i = 1:numberOfEra
%        figure(i);
%        drawVariability(eraDistanceScore(i,:,:) , onlyEraDistanceScore(i,:,:), agesDescription{i},'region index');
%     end
end

function drawVariability(eraDistanceScore , onlyEraDistanceScore, regionName,xtitle)
%         subplot(1,2,1);
% %         plot( squeeze(mean(eraDistanceScore,2)) );
% %         errorbar(squeeze(mean(eraDistanceScore,2)), squeeze(std(eraDistanceScore,1,2)) ,'.');
%         plotp(eraDistanceScore', 1:size(eraDistanceScore,1), 'r-', false);
% 
%         title(sprintf('%s scores without genes', regionName),'fontsize',16);
%         xlabel(xtitle,'fontsize',20);
%         ylabel('variability score (between/within)','fontsize',20);
%         ylim([5 12]);
%         subplot(1,2,2);
%         plot( squeeze(mean(onlyEraDistanceScore,2)) );
%         errorbar(squeeze(mean(onlyEraDistanceScore,2)), squeeze(std(onlyEraDistanceScore,1,2)) ,'.');
        plotp(onlyEraDistanceScore');
        title(sprintf('%s scores using only the  genes', regionName),'fontsize',16);
        xlabel(xtitle,'fontsize',20);
        ylabel('variability score (between/within)','fontsize',20);
        ylim([5 12]);
end

