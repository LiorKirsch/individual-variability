function [expression, gross_region_vec, gene_symbols, samples2subjects, gross_structures_names] = load_expression_and_regions(dataname, parms)

 switch dataname
      case 'human6', 
        [expression, gross_region_vec, gene_symbols, samples2subjects, gross_structures_names] = load_human6(parms);
      case 'kang', 
        [expression, gross_region_vec, gene_symbols, samples2subjects, gross_structures_names] = load_kang(parms);
      case 'mouse', 
        [expression, gross_region_vec, gene_symbols, samples2subjects, gross_structures_names] = load_mouse(parms);
    end

end

% ===========================================================
function  [expression, gross_region_vec, gene_symbols,samples2subjects,gross_structures_names] = load_human6( parms)
% 
    
    persistent local_experimentsDataMatrix
    persistent local_experimentsLocationMatrix    
    persistent local_humanOntology
    persistent local_gene_symbols
    persistent local_samples2subjects
    persistent local_grossStructures

    if isempty(local_experimentsDataMatrix)
        fprintf('human6: loading data from disk\n');
        % Load human data 
        data_dirname = fullfile('/','home','lab', 'lior', 'Projects', 'individual variability');
        filename = 'easyFormatHumanData.mat';
        fullname = fullfile(data_dirname, filename);
        load(fullname, 'experimentsLocationMatrix', 'experimentsDataMatrix','selectedProbesData','experimentsSubjectMatrixLogical'); %#ok
        
        % Load the tree
        filename = 'humanOntologyObject.mat';
        load(fullfile(data_dirname, filename), 'humanOntology'); %#ok
        
        grossStructures = {'Frontal Lobe';'Cingulate gyrus';'hippocampal formation';'parahippocampal gyrus';...
                           'Occipital Lobe';'Parietal Lobe';'Temporal Lobe';'Amygdala';'Basal Forebrain';...
                           'Striatum';'Hypothalamus';'Dorsal Thalamus';'Mesencephalon';'Cerebellar Cortex';...
                           'Pontine Tegmentum';'Myelencephalon'};
                       
        local_humanOntology = humanOntology;
        local_experimentsDataMatrix = experimentsDataMatrix;
        local_experimentsLocationMatrix = experimentsLocationMatrix;
        local_gene_symbols = selectedProbesData.gene_symbols;
        local_samples2subjects = experimentsSubjectMatrixLogical;
        local_grossStructures = grossStructures;
    end
    
    gross_region_vec = fine_to_gross(local_experimentsLocationMatrix, local_humanOntology,local_grossStructures);
    relevant_samples = (gross_region_vec>0);    
    expression = local_experimentsDataMatrix(relevant_samples,:);
    gross_region_vec = gross_region_vec(relevant_samples);
    gene_symbols = local_gene_symbols;
    samples2subjects = local_samples2subjects(relevant_samples,:);
    gross_structures_names = local_grossStructures;
end

% ===========================================================
function  [expression, gross_region_vec, gene_symbols, samples2subjects,gross_structures_names] = load_mouse(parms)
% 
    persistent local_expression;
    persistent local_gross_region_vec    
    persistent local_gene_symbols;
    persistent local_grossStructures
    
    if isempty(local_expression)
        fprintf('mouse: loading data from disk\n');
        % Load human data 
        data_dirname = fullfile('/', 'home','lab', 'gal', 'Projects', ...
                                'Limor', 'Data');
        filename = 'ABA_expression_data.mat';
        fullname = fullfile(data_dirname, filename);
        x = load(fullname);
        
        % Focuse on gross structures
        grossStructures = {'Cerebral cortex'; 'Olfactory areas'; ...
                           'Hippocampal region'; 'Retrohippocampal region'; ...
                           'Striatum'; 'Pallidum'; ...
                           'Hypothalamus'; 'Thalamus'; 'Midbrain'; ...
                           'Cerebellum'; 'Pons'; 'Medulla'};
        num_gross = numel(grossStructures);
        
        inds = cellfind(x.name_of_regions, grossStructures);
        if any(inds==0), error('Name or region not found');end
        expression = x.expression_energy(inds,:);
        gross_region_vec = (1:num_gross);
        
        local_expression = expression;
        local_gross_region_vec = gross_region_vec(:);
        local_gene_symbols = x.name_of_genes;
        local_grossStructures = grossStructures;
    end
    expression = local_expression;
    gross_region_vec = local_gross_region_vec(:);
    gene_symbols = local_gene_symbols;
    samples2subjects = ones(size(gross_region_vec)); % the mouse does not have different subjects;
    gross_structures_names = local_grossStructures;
end



% ===========================================================
function  [expression, gross_region_vec,gene_symbols,gross_structures_names] = load_kang(parms)
% 
    persistent local_expression;
    persistent local_gross_region_vec    
    persistent local_gene_symbols;
    persistent local_grossStructures

    if isempty(local_expression)
        fprintf('kang: loading data from disk\n');
        % Load human data 
        data_dirname = fullfile('/', 'cortex','data', 'microarray', 'human', ...
                                'Kang2011');
        filename = 'kang_individual_samples.mat';
        fullname = fullfile(data_dirname, filename);
        x = load(fullname);
        
        % Focuse on gross structures
        %
        % We need to map between 4 layes: 
        % samples->regions->gross->sorted-strurctures
 
        
        % Map from gross to sorted
        grossStructures = { 'Frontal Lobe','Hippocampus',  ...
                            'Occipital Lobe', 'Parietal Lobe', ...
                            'Temporal Lobe', 'Amygdala', 'Striatum', ...
                            'Thalamus', 'Cerebellum'};        
        sorted_inds = cellfind(x.grossRegionNames, grossStructures);
        assert( all(0 < sorted_inds ) ,'Name of region not found');
        
        % Map from samples to gross: 
        grossRegionSamples = double(x.samples2regions) * double(x.regionOntology');
        samples2grossID = grossRegionSamples * sorted_inds';
        is_sample_good_regions = (samples2grossID>0);
 
      
        % Focuse on adult periods        
        x.periods = x.samples2periods * (1:size(x.samples2periods, 2))';
        adults = {'12-20y', '20-40y', '40-60y', '60+y'};        
        is_sample_good_ages = ismember(x.period_names(x.periods), adults);
        
        % Focus on adults and on gross regions
        relevent_samples = is_sample_good_regions & is_sample_good_ages';
        expression = x.data(:, relevent_samples)';
        gross_region_vec =  samples2grossID( relevent_samples );
        
        
        local_expression = expression;
        local_gross_region_vec = gross_region_vec(:);
        local_gene_symbols = x.gene_names;
        local_grossStructures = grossStructures;

    end
    expression = local_expression;
    gross_region_vec = local_gross_region_vec(:);
    gene_symbols = local_gene_symbols;
    gross_structures_names = local_grossStructures;

end
