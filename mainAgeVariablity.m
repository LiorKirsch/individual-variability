load('easyFormatHumanData_withAges.mat');


load('subset.mat');
load('humanOntologyObject.mat');
regionColors = humanOntology.structureColors(subsetNodesIndex,:);
regionNames = humanOntology.structureLabels(subsetNodesIndex,4);
experimentRegion = getRegionOfGrossStructures(subsetNodesIndex, experimentsLocationMatrix, humanOntology);

numberOfRegions = size(experimentRegion,2);
numberOfPeople = size(experimentsSubjectMatrixLogical,2);
 
numberOfSamples = size(experimentsSubjectMatrixLogical,1);
numberOfRegions = size(experimentRegion,2);
 
assert(  numberOfSamples == size(experimentsDataMatrix,1) );

allSamples = true(size(experimentRegion,1),1);
[categoriesScoresAllSamples, distanceBetweenGroupsAllSamples] = calcRegionBetweenSubjectDist(experimentsSubjectMatrixLogical, experimentsDataMatrix, allSamples);
numberOfEra = size(agesDescription,1);

eraDistanceScoreAllSamples = nan(numberOfEra,numberOfPeople);
onlyEraDistanceScoreAllSamples = nan(numberOfEra,numberOfPeople);
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

[categoriesScores, distanceBetweenGroups] = calcRegionBetweenSubjectDist(experimentsSubjectMatrixLogical, experimentsDataMatrix, experimentRegion);

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


figure(12121);

subplot(1,2,1);
%plot( squeeze(mean(eraDistanceScoreAllSamples,2)) );
errorbar(squeeze(mean(eraDistanceScoreAllSamples,2)), squeeze(std(eraDistanceScoreAllSamples,1,2)) ,'.');
title('all regions - scores without age genes');
subplot(1,2,2);
% plot( squeeze(mean(onlyEraDistanceScoreAllSamples,2)) );
errorbar(squeeze(mean(onlyEraDistanceScoreAllSamples,2)), squeeze(std(eraDistanceScoreAllSamples,1,2)) ,'.');
title('all regions - scores with only the age genes');
% xticklabel_rotate(1:length(agesDescription),40,agesDescription, 'fontSize',14) ;

for m = 1:numberOfRegions
   figure(m);
   
    subplot(1,2,1);
    %plot( squeeze(mean(eraDistanceScoreAllSamples,2)) );
    errorbar(squeeze(mean(eraDistanceScore(:,:,m),2)), squeeze(std(eraDistanceScore(:,:,m),1,2)) ,'.');
    title(sprintf('%s scores without age genes', regionNames{m}));
    subplot(1,2,2);
    % plot( squeeze(mean(onlyEraDistanceScoreAllSamples,2)) );
    errorbar(squeeze(mean(onlyEraDistanceScore(:,:,m),2)), squeeze(std(eraDistanceScore(:,:,m),1,2)) ,'.');
    title(sprintf('%s - scores with only the age genes', regionNames{m}));
   
end
