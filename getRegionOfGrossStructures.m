function experimentInGrossRegion = getRegionOfGrossStructures(grossStructureIndicesInOntlogy, experimentLocationInOntology, ontologyObject)
% experimentInGrossRegion = getRegionOfGrossStructures(grossStructureIndicesInOntlogy, experimentLocationInOntology, ontologyObject)
% generate a logical matrix which experiment belong to which of the gross
% area provided in grossStructureIndicesInOntlogy
% grossStructureIndicesInOntlogy - logical array, hold the gross area indices
% experimentLocationInOntology - for each experiment it hold the location in the ontology
% ontologyObject - the ontology object
    
    numOfRegionInOntology = length(grossStructureIndicesInOntlogy);
    numOfGrossRegions = sum(grossStructureIndicesInOntlogy);
    numOfExperiments = size(experimentLocationInOntology,1);
    
    assert( numOfRegionInOntology == size(experimentLocationInOntology,2) );
        
    [allChilds,~] = ontologyObject.allChildNodes();
    grossRegionChilds = allChilds(grossStructureIndicesInOntlogy,:);
    
    experimentInGrossRegion =  false(numOfExperiments, numOfGrossRegions);
    for i=1:numOfGrossRegions % loop over the high category areas
        currentGrossAreaChilds = grossRegionChilds(i,:)';
        releventExperiments = double(experimentLocationInOntology) * double(currentGrossAreaChilds);
        experimentInGrossRegion(:,i) = logical(releventExperiments);
    end
end