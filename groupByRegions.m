function [experimentsDataMatrix, locationMatrix, experimentsSubjectMatrixLogical, experimentRegion, mni_xyz, mri_voxel_xyz]= groupByRegions(experimentsDataMatrix, locationMatrix, experimentsSubjectMatrixLogical,ontologyObject, subsetNodesIndex, mni_xyz, mri_voxel_xyz)

    assert(size(experimentsDataMatrix,1) == size(locationMatrix,1));
    assert(size(locationMatrix,2) == size(ontologyObject.dependencyMatrix,1) );
    
    experimentRegion = getRegionOfGrossStructures(subsetNodesIndex, locationMatrix, ontologyObject);

    sampleWithNoLocation = ~any(experimentRegion,2);
    experimentsDataMatrix = experimentsDataMatrix(~sampleWithNoLocation,:);
    locationMatrix = locationMatrix(~sampleWithNoLocation,:);
    experimentsSubjectMatrixLogical = experimentsSubjectMatrixLogical(~sampleWithNoLocation,:);
    experimentRegion = experimentRegion(~sampleWithNoLocation,:);
    mni_xyz = mni_xyz(~sampleWithNoLocation,:);
    mri_voxel_xyz = mri_voxel_xyz(~sampleWithNoLocation,:);
end