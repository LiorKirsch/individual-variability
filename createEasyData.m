function createEasyData()

load('expressionRawDataWithOntology.mat', 'expressionDataSelected',  'humanOntology','selectedProbesData','locationMatrix','structureData');
load('humanOntologyObject.mat');
numberOfPeople = 6;

[experimentsDataMatrix, experimentsLocationMatrix, experimentsSubjectMatrix, probeCode, mni_xyz, mri_voxel_xyz] = joinDataToOneBigMatrix(expressionDataSelected, locationMatrix, structureData);

not_a_probe = ~strcmp(selectedProbesData.chromosome,'');

experimentsDataMatrix = experimentsDataMatrix(:, not_a_probe);

selectedProbesData.probe_ids = selectedProbesData.probe_ids(not_a_probe);
selectedProbesData.probe_names = selectedProbesData.probe_names(not_a_probe);
selectedProbesData.gene_ids = selectedProbesData.gene_ids(not_a_probe);
selectedProbesData.gene_names = selectedProbesData.gene_names(not_a_probe);
selectedProbesData.gene_symbols = selectedProbesData.gene_symbols(not_a_probe);
selectedProbesData.entrez_ids = selectedProbesData.entrez_ids(not_a_probe);
selectedProbesData.chromosome = selectedProbesData.chromosome(not_a_probe);
selectedProbesData.bestProbeForGene = selectedProbesData.bestProbeForGene(not_a_probe);

    
experimentsSubjectMatrixLogical = false(size(experimentsSubjectMatrix,1), numberOfPeople);
for i =1:numberOfPeople
    experimentsSubjectMatrixLogical(:,i) = experimentsSubjectMatrix == i;
end

save('easyFormatHumanData.mat','experimentsDataMatrix', 'experimentsLocationMatrix', 'experimentsSubjectMatrixLogical', 'selectedProbesData', 'mni_xyz', 'mri_voxel_xyz');
end

