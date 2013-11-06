addpath('/home/lab/lior/Projects/datasets/Human allen brain/@ontology');


% Load data 
dirname = fullfile('/', 'home', 'lab', 'lior', 'Projects', ...
		   'human variability');
filename = 'dataMatrix24.mat'
x = load(fullfile(dirname, filename));

[num_genes, num_regions, num_subjects] = size(x.dataMatrix);
% datamatrix is the mean over all small regions, for each
% super-region.  after first averaging over samples within a small
% region.

% Compute variability for each region

