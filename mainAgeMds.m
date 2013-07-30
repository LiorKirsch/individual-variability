
load('easyFormatHumanData_withAges.mat');
% load('humanOntologyObject.mat');
% load('subset.mat');

% remove genes that have more (or less) then one age
validGenes = sum(ages,2) == 1 ;
experimentsDataMatrix = experimentsDataMatrix(:, validGenes);
ages = logical(ages(validGenes,:));

ageColors = createColorMap(length(agesDescription));
plotMds(experimentsDataMatrix',ages, NaN, agesDescription,ageColors);

