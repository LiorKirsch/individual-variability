function [experimentsDataMatrix, experimentsLocationMatrix, experimentsSubjectMatrix, probeCode, mni_xyz, mri_voxel_xyz] = joinDataToOneBigMatrix(experimentData, locationMatrix, structureData)
    numberOfProbes = size(experimentData{1}.expressionMatrix, 1);
    numberOfSamples = nan(length(experimentData),1);
    numberOfSubjects = length(experimentData);
    numberOfAreas = size(locationMatrix{1}, 2);
    
    for i=1:numberOfSubjects
        numberOfSamples(i) = size(experimentData{i}.expressionMatrix, 2);
    end
    totalNumberOfSamples = sum(numberOfSamples);
    
    experimentsDataMatrix = nan(totalNumberOfSamples, numberOfProbes);
    experimentsLocationMatrix = nan(totalNumberOfSamples, numberOfAreas);
    experimentsSubjectMatrix = nan(totalNumberOfSamples,1);
    
    mni_xyz = nan(totalNumberOfSamples, 3);
    mri_voxel_xyz = nan(totalNumberOfSamples, 3);

    lastIndex = 1;
    probeCode = experimentData{1}.probeCode;
    for i = 1:numberOfSubjects
       assert(all(probeCode == experimentData{i}.probeCode));
       indicesOfIntrests = lastIndex: lastIndex + numberOfSamples(i) - 1;
       experimentsDataMatrix(indicesOfIntrests,:) = experimentData{i}.expressionMatrix';
       experimentsLocationMatrix( indicesOfIntrests , :) = locationMatrix{i};
       experimentsSubjectMatrix ( indicesOfIntrests) =  i ;
       mni_xyz(indicesOfIntrests,:) = structureData{i}.mni_xyz;
       mri_voxel_xyz(indicesOfIntrests,:) = structureData{i}.mri_voxel_xyz;
       lastIndex = lastIndex + numberOfSamples(i);
    end
        
end