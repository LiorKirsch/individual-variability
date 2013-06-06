load('expressionRawDataWithOntology.mat', 'expressionDataSelected',  'humanOntology','selectedProbesData','locationMatrix','structureData');

load('humanOntologyObject.mat');
load('subset.mat');

[aboveK ,numOfSamples] = checkForMoreThenKSamples(locationMatrix, humanOntology,10);
subsetNodesIndex = subsetNodesIndex & aboveK;
%subsetNodesIndex = true(size(aboveK));
regionColors = humanOntology.structureColors(subsetNodesIndex,:);
regionNames = humanOntology.structureLabels(subsetNodesIndex,4);



[experimentsDataMatrix, experimentsLocationMatrix, experimentsSubjectMatrix, probeCode, mni_xyz, mri_voxel_xyz] = joinDataToOneBigMatrix(expressionDataSelected, locationMatrix,structureData);


[experimentsDataMatrix, experimentsLocationMatrix, experimentsSubjectMatrix, probeCode, experimentRegion] = calcPCAinStructures(expressionDataSelected, locationMatrix, humanOntology,subsetNodesIndex);
[coeff,score,latent] = pca(experimentsDataMatrix);
save('pcaData.mat','experimentsDataMatrix', 'experimentsLocationMatrix', 'experimentsSubjectMatrix', 'probeCode', 'experimentRegion', 'coeff','score', 'latent','regionNames','regionColors','mni_xyz', 'mri_voxel_xyz');
%save('pcaDataAllRegions.mat','experimentsDataMatrix', 'experimentsLocationMatrix', 'experimentsSubjectMatrix', 'probeCode', 'experimentRegion', 'coeff','score', 'latent','regionNames','regionColors');


% % some checks that everything is in the same order
% locationString = cell(size(experimentsLocationMatrix,1),1);
% for i =1:size(experimentsLocationMatrix,1)
%     expIndices = logical(experimentsLocationMatrix(i,:));
%     locationString{i} = humanOntology.structureLabels( expIndices,4);
%     
% end