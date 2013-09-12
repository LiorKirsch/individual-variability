close all
% load('easyFormatHumanData_withAges.mat');
experimentsDataMatrixNormalized = normalizeMatrixByDim(experimentsDataMatrix, 1);
numberoOfSamples = size(experimentsDataMatrix,1);
numberoOfGenes = size(experimentsDataMatrix,2);
numberOfEra = size(agesDescription,1);



load('subset.mat');
load('humanOntologyObject.mat');
% subsetNodesBooleanIndex =  [5         122         145         155         200         248         303         752         925        1114 ];
% subsetNodesIndex(:) = false;
% subsetNodesIndex(subsetNodesBooleanIndex) = true;

regionNames = humanOntology.structureLabels(subsetNodesIndex,4);
experimentRegion = getRegionOfGrossStructures(subsetNodesIndex, experimentsLocationMatrix, humanOntology);

regionColors = humanOntology.getColorByRegionName(regionNames);

numberOfRegions = size(experimentRegion,2);


genesWithMoreThenOneAge = sum(ages,2) > 1;
agesNumbers = agesWithoutTheseGenes * (1:size(ages,2) )';
agesNumbers(genesWithMoreThenOneAge) = 0;

TAI_scores = nan(numberOfRegions,1);
TAI_scores_std = nan(numberOfRegions,1);

samples_TAI_scores = (experimentsDataMatrix * agesNumbers ) ./ sum(experimentsDataMatrix,2);

for i = 1:numberOfRegions
    
   regionName =  regionNames{i};
   experimentsInRegion = experimentRegion(:,i);
   
   TAI_scores(i) = mean(samples_TAI_scores(experimentsInRegion) );
   TAI_scores_std(i) = std(samples_TAI_scores(experimentsInRegion) );
   

end

drawBars(TAI_scores,TAI_scores_std,regionNames,regionColors,'TAI');









TAI_scores = nan(numberOfRegions,numberOfEra);
TAI_scores_std = nan(numberOfRegions,1);

samples_TAI_scores = (experimentsDataMatrix * agesNumbers ) ./ sum(experimentsDataMatrix,2);
normSamples_scores = diag( 1./ sum(experimentsDataMatrix,2) ) * experimentsDataMatrix  ;


for i = 1:numberOfRegions

   regionName =  regionNames{i};
   experimentsInRegion = experimentRegion(:,i);

   TAI_scores_std(i) = std(samples_TAI_scores(experimentsInRegion) );
   for j = 1:numberOfEra
        eraGenes = (agesNumbers == j);
        scores = sum(normSamples_scores(experimentsInRegion,eraGenes) * j , 2);
        TAI_scores(i,j) = mean(scores );
   end
end
 bar(TAI_scores, 'stack');