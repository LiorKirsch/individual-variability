function mainCalcSparsePca()


    load('easyFormatHumanData.mat');
    load('humanOntologyObject.mat');
    load('subset.mat');

    [aboveK ,numOfSamples] = checkForMoreThenKSamples(experimentsLocationMatrix, experimentsSubjectMatrixLogical, humanOntology,5);
    subsetNodesIndex = subsetNodesIndex & aboveK;
    %subsetNodesIndex = true(size(aboveK));
    regionColors = humanOntology.structureColors(subsetNodesIndex,:);
    regionNames = humanOntology.structureLabels(subsetNodesIndex,4);

    [experimentsDataMatrix, experimentsLocationMatrix, experimentsSubjectMatrix, experimentRegion, mni_xyz, mri_voxel_xyz] = groupByRegions(experimentsDataMatrix, experimentsLocationMatrix, experimentsSubjectMatrixLogical,humanOntology,subsetNodesIndex, mni_xyz, mri_voxel_xyz);
    numberOfSamples = size(experimentsDataMatrix,1);

    m=4;                                % Number of components:
    gamma=0.5*ones(1,m);                % sparsity weight factors -one for each component - 
                                        % in relative value with respect to the theoretical upper bound
    
    experimentsDataMatrix = experimentsDataMatrix-repmat((mean(experimentsDataMatrix,1)),numberOfSamples,1);        % Centering of the data

    coeff = GPower(experimentsDataMatrix,gamma,m,'l1',0); 
    score = experimentsDataMatrix * coeff;
    latent = 1;

    save('pcaSparseData.mat','experimentsDataMatrix', 'experimentsLocationMatrix', 'experimentsSubjectMatrix', 'selectedProbesData', 'experimentRegion', 'coeff','score', 'latent','regionNames','regionColors','mni_xyz', 'mri_voxel_xyz');

end