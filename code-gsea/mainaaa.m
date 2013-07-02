% load('expressionRawDataWithOntology.mat', 'expressionDataSelected',  'humanOntology','selectedProbesData','locationMatrix');
% load('subset.mat')
% [anovaMatrix, samplesInRegion] = calcAnovaInStructures(expressionDataSelected, locationMatrix, humanOntology,areas24indices);
% regionLabels = humanOntology.structureLabels(areas24indices , 4);
% save('anovaResults.mat','anovaMatrix', 'selectedProbesData', 'regionLabels','samplesInRegion');


load('anovaResults.mat','anovaMatrix', 'selectedProbesData', 'regionLabels','samplesInRegion');

for i = 1:length(regionLabels)
    outputFile = ['listOfGenes/' , regionLabels{i} ,'.txt'];
    [scores, indices] = sort(anovaMatrix(:,i));
    sortedGeneNames = selectedProbesData.gene_symbols(indices);
    fileID = fopen(outputFile,'w');
    for j=1:length(sortedGeneNames)
        fprintf(fileID,'%s\r\n',sortedGeneNames{j});
    end
end

%meanOfRegions = mean(log(anovaMatrix),1)';
%regionLabels = humanOntology.structureLabels(areas24indices, 4);
%regionLabelsSorted = regionLabels(indices);



% fileID = fopen('regionScores.csv','w');
% fprintf(fileID,'mean(log(anova)),region,#samples1,#samples2,#samples3,#samples4,#samples5,#samples6\r\n');
% for j=1:length(regionLabels)
%     fprintf(fileID,'%g,%s,%d, %d, %d, %d, %d, %d \r\n',meanOfRegions(j,:), regionLabels{j}, samplesInArea(j,:));
% end 

%genesInCategory = indicator of gene in category

geneGo = load('humanGene2GoMatrix.mat', 'go_genes_mat', 'cat_ids', 'geneNames', 'aspects');
