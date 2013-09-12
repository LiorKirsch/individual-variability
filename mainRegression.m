
load('easyFormatHumanData.mat');
load('humanOntologyObject.mat');
load('subset.mat');

[aboveK ,numOfSamples] = checkForMoreThenKSamples(experimentsLocationMatrix, experimentsSubjectMatrixLogical, humanOntology,5);
subsetNodesIndex = subsetNodesIndex & aboveK;
%subsetNodesIndex = true(size(aboveK));
regionColors = humanOntology.structureColors(subsetNodesIndex,:);
regionNames = humanOntology.structureLabels(subsetNodesIndex,4);

[experimentsDataMatrix, experimentsLocationMatrix, experimentsSubjectMatrix, experimentRegion, mni_xyz, mri_voxel_xyz] = groupByRegions(experimentsDataMatrix, experimentsLocationMatrix, experimentsSubjectMatrixLogical,humanOntology,subsetNodesIndex, mni_xyz, mri_voxel_xyz);



experimentRegionIndecies = experimentRegion * (1:size(experimentRegion,2))';
% [linearRegressionWeights, lassoWeights, ordinalWeights] = findWeightsOfFit(experimentsDataMatrix, experimentRegionIndecies);

