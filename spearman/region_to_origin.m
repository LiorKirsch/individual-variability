function gross_mat = region_to_origin(region_mat, ontologyObject)
% 
% map fine regions (1:1534) to coarse regions (1:10,NaN)
%
    
    [num_experiments, num_fine] = size(regions_mat);
    
    gross_regions_inds = [1 ,5 ,6, 10, 1534]
    num_gross = length(gross_regions_inds);    
    gross_regions = sparse(logical(1, num_fine));
    gross_regions(gross_regions_inds) = true;
    
        
    allChilds = ontologyObject.allChildNodes();
    grossRegionChilds = allChilds(gross_regions, :);
    
    experimentInGrossRegion =  false(numOfExperiments, numOfGrossRegions);
    for i=1:numOfGrossRegions % loop over the high category areas
        currentGrossAreaChilds = grossRegionChilds(i,:)';
        releventExperiments = double(experimentLocationInOntology) * double(currentGrossAreaChilds);
        experimentInGrossRegion(:,i) = logical(releventExperiments);
    end
end
~
~

    
    
end