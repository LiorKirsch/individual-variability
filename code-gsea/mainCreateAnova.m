clear; 
load('expressionRawDataWithOntology.mat', 'expressionDataSelected',  'humanOntology','selectedProbesData','locationMatrix');
load('subset.mat')
[anovaMatrix, samplesInRegion] = calcAnovaInStructures(expressionDataSelected, locationMatrix, humanOntology,areas24indices);
regionLabels = humanOntology.structureLabels(areas24indices , 4);
save('anovaResults.mat','anovaMatrix', 'selectedProbesData', 'regionLabels','samplesInRegion');
