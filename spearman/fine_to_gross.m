function gross_region_vec = fine_to_gross(fine_region_mat, ontologyObject,grossStructures)
% 
% map fine regions (1:1534) to coarse regions (1:10,NaN)
% fine_region_mat - list of indcies of the finer region
% ontlogyObject - contains the mapping between the regions
% grossStructures - are the names of the gross regions
% 
% Output: 
%   returns a vector mapping eavh of the 3702 regions to a number in (1:numel(grossStructures)) 
    
    [num_experiments, num_fine] = size(fine_region_mat);
    all_structures_names = ontologyObject.structureLabels(:,4);
    
    [are_members, gross_regions_inds] = ismember(grossStructures, all_structures_names);
    assert(all(are_members),'not all grossStructures appears in the ontology');

    num_gross = length(gross_regions_inds);    
%     gross_regions = logical(sparse(1, num_fine));
%     gross_regions(gross_regions_inds) = true;
    
        
    allChilds = ontologyObject.allChildNodes();
    grossRegionChilds = allChilds(gross_regions_inds, :);
    
    gross_region_mat = false(num_experiments, num_gross);
    for i = 1:num_gross % loop over the gross regions
        currentGrossAreaChilds = grossRegionChilds(i,:)';
        releventExperiments = double(fine_region_mat) ...
            * double(currentGrossAreaChilds);
        gross_region_mat(:,i) = logical(releventExperiments);
    end
    gross_region_vec = gross_region_mat * (1:num_gross)';
end
