close all
% load('easyFormatHumanData_withAges.mat');
% experimentsDataMatrixNormalized = experimentsDataMatrix - repmat(mean(experimentsDataMatrix,1), [size(experimentsDataMatrix,1), 1] ) * diag(1 ./ std(experimentsDataMatrix,1));


numberOfEra = size(agesDescription,1);
expressionOfEachSamples = nan(size(experimentsDataMatrix,1), numberOfEra);
for i = 1:numberOfEra
   
   genesInEra = logical(ages(:,i));
    
   currentDataMatrix = experimentsDataMatrix(:,genesInEra);
   expressionOfEachSamples(:,i) = mean(currentDataMatrix,2);
end


expressionOfEachSamples = diag( 1 ./ mean(expressionOfEachSamples,2) ) * expressionOfEachSamples;

% for i = 1:numberOfEra
%    genesInEra = logical(ages(:,i));
%    eraName = agesDescription{i}; 
%    figure(i);
%    scatter3(mri_voxel_xyz(:,1), mri_voxel_xyz(:,2) ,mri_voxel_xyz(:,3) ,5 ,expressionOfEachSamples(:,i) );
%    
%    titleString = sprintf('%s - (%d genes) ' , eraName, sum(genesInEra));
%    title(titleString); 
% end




load('subset.mat');
load('humanOntologyObject.mat');
regionColors = humanOntology.structureColors(subsetNodesIndex,:);
regionNames = humanOntology.structureLabels(subsetNodesIndex,4);
experimentRegion = getRegionOfGrossStructures(subsetNodesIndex, experimentsLocationMatrix, humanOntology);

numberOfRegions = size(experimentRegion,2);

% for i = 1:numberOfRegions
%     
%    regionName =  regionNames{i};
%    experimentsInRegion = experimentRegion(:,i);
%    
%    meanExpressionInEras = nan(numberOfEra,1);
%    stdExpressionInEras = nan(numberOfEra,1);
%    for j = 1:numberOfEra
%        genesInEra = logical(ages(:,j));
% 
%        currentDataMatrix = experimentsDataMatrix(experimentsInRegion,genesInEra);
%        expressionOfSamplesInRegion = mean(currentDataMatrix,2);
%        meanExpressionInEras(j) = mean(expressionOfSamplesInRegion);
%        stdExpressionInEras(j) = std(expressionOfSamplesInRegion);
%    end
%    
%    figure(i*10);
%    hold on;
% %    bar(meanExpressionInEras);
%    errorbar(meanExpressionInEras, stdExpressionInEras, '.');
%    ylim([3.5 7.5]);
%    ylabel('mean expression', 'fontSize',20);
%    xlabel('gene age (era)', 'fontSize',20);
%    titleString = sprintf('%s - (%d samples) ' , regionName, sum(experimentsInRegion));
%    title(titleString, 'fontSize',20);
% end

%%
for j = 1:numberOfEra

   eraName = agesDescription{j}; 
   
   genesInEra = logical(ages(:,j));
   
   
   
   meanExpressionInEras = nan(numberOfEra,1);
   stdExpressionInEras = nan(numberOfEra,1);
   
   for i = 1:numberOfRegions
       
       experimentsInRegion = experimentRegion(:,i);
       currentDataMatrix = experimentsDataMatrix(experimentsInRegion,genesInEra);
%       currentDataMatrix = experimentsDataMatrixNormalized(experimentsInRegion,genesInEra);
       expressionOfSamplesInRegion = mean(currentDataMatrix,2);
       meanExpressionInEras(i) = mean(expressionOfSamplesInRegion);
       stdExpressionInEras(i) = std(expressionOfSamplesInRegion);
   end
   
   figure(j*100);
   hold on;
%    bar(meanExpressionInEras);
   errorbar(meanExpressionInEras, stdExpressionInEras, '.');
   ylabel('mean expression', 'fontSize',20);
   xlabel('brain region', 'fontSize',20);
   ylim([3.5 7.5]);
   titleString = sprintf('%s - (%d genes) ' , eraName, sum(genesInEra));
   title(titleString, 'fontSize',20);
   xticklabel_rotate(1:length(regionNames),40,regionNames, 'fontSize',14) ;
end


%%

aucScores = nan(numberOfEra, numberOfRegions);
for j = 1:numberOfEra

   eraName = agesDescription{j}; 
   
   genesInEra = logical(ages(:,j));
   
   for i = 1:numberOfRegions
       
       experimentsInRegion = experimentRegion(:,i);
       currentDataMatrix = experimentsDataMatrix(experimentsInRegion,:);
%       currentDataMatrix = experimentsDataMatrixNormalized(experimentsInRegion,:);
       meanExpressionInArea = mean(currentDataMatrix,1)';
       
       [aucScores(j,i),~,~] = myauc(genesInEra, meanExpressionInArea,true);
       
   end
   
  
end

imagesc(aucScores);
colormap(jet); colorbar; ylabel('ages'); xlabel('region index');