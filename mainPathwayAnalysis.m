function mainPathwayAnalysis()

    load('gene_pathways.mat');
    load('easyFormatHumanData_withAges.mat');
    load('subset.mat');
    load('humanOntologyObject.mat');
    regionSampleTreshold = 50;

    [transmitters_genes, transmitters] = getTransmitorsGenes(path_genes_mat, selectedProbesData);
    
    numberOfTransmitters = length( transmitters );
    numberOfPeople = size(experimentsSubjectMatrixLogical,2);
    numberOfSamples = size(experimentsSubjectMatrixLogical,1);

    assert(  numberOfSamples == size(experimentsDataMatrix,1) );
    

    regionNames = humanOntology.structureLabels(subsetNodesIndex,4);
    experimentRegion = getRegionOfGrossStructures(subsetNodesIndex, experimentsLocationMatrix, humanOntology);

    experimentsInRegion = sum(experimentRegion,1);
    aboveTreshold = (experimentsInRegion > regionSampleTreshold);
    experimentRegion = experimentRegion(:,aboveTreshold);
    regionNames = regionNames(aboveTreshold);
    regionColors = humanOntology.getColorByRegionName(regionNames);
    experimentsInRegion = experimentsInRegion(aboveTreshold);
    
    
%     experimentsDataMatrixNormalized = normalizeMatrixByDim(experimentsDataMatrix, 1);
%     [onlySubsetDistanceScore, withoutSubsetDistanceScore] = calcClusteringScoreWithAndWithoutSubset(experimentsSubjectMatrixLogical, experimentsDataMatrixNormalized, allSamples, transmitters_genes);
    [onlySubsetDistanceScore, withoutSubsetDistanceScore] = calcClusteringScoreWithAndWithoutSubsets(experimentsSubjectMatrixLogical, experimentsDataMatrix, experimentRegion, transmitters_genes);

    colors = {'r','b','g','k','y'};
    f=figure;hold on;
    for i = 1:numberOfTransmitters
       data = squeeze(onlySubsetDistanceScore(i,:,:) );
       plot( mean(data,1) ,colors{i},'LineWidth',1.5);
       %drawVariability(data, colors{i}, 'only genes','region index',f);
    end
    xticklabel_rotate(1:length(regionNames),40,regionNames, 'fontSize',14) ;
    legend(transmitters);
    hold off;
    
    data2 = squeeze(mean(onlySubsetDistanceScore,2));
    bar(data2');
    xticklabel_rotate(1:length(regionNames),40,regionNames, 'fontSize',14) ;
    legend(transmitters);
end

function drawVariability(data, color,theTitle,xtitle,f)
        plotp(1:size(data,2),data,color,f);
        title( theTitle ,'fontsize',16);
        xlabel(xtitle,'fontsize',20);
        ylabel('variability score (between/within)','fontsize',20);
end

function [onlySubsetDistanceScore, withoutSubsetDistanceScore] = calcClusteringScoreWithAndWithoutSubsets(experimentsSubjectMatrixLogical, experimentsDataMatrix, samplesSubset, geneSubsets)

    
    numberOfPeople = size(experimentsSubjectMatrixLogical,2);
    numberOfGeneSubsets = size(geneSubsets,1);
    numberOfGenes = size(experimentsDataMatrix,2);
    numberOfSampleSubsets = size(samplesSubset,2);

    assert( numberOfGenes == size(geneSubsets,2));
    
  
    
    [categoriesScores, distanceBetweenGroups] = calcRegionBetweenSubjectDist(experimentsSubjectMatrixLogical, experimentsDataMatrix, samplesSubset);
    
    withoutSubsetDistanceScore = nan(numberOfGeneSubsets,numberOfPeople,numberOfSampleSubsets);
    onlySubsetDistanceScore = nan(numberOfGeneSubsets,numberOfPeople,numberOfSampleSubsets);
    
    for i = 1:numberOfGeneSubsets

        genesInSubset = geneSubsets(i,:);
        currentDistanceBetweenGroups = distanceBetweenGroups(:,:,:,~genesInSubset);
        currentDistanceBetweenGroupsOnlyEra = distanceBetweenGroups(:,:,:,genesInSubset);

        for m =1:numberOfSampleSubsets
            currentAreaDistanceBetweenGroups = squeeze(currentDistanceBetweenGroups(m,:,:,:));
            withoutSubsetDistanceScore(i,:,m) = clusteringScore( currentAreaDistanceBetweenGroups );

            currentAreaOnlyCatDistanceBetweenGroups = squeeze(currentDistanceBetweenGroupsOnlyEra(m,:,:,:));
            onlySubsetDistanceScore(i,:,m) = clusteringScore( currentAreaOnlyCatDistanceBetweenGroups );
        end
        printPercentCounter(i, numberOfGeneSubsets );
    end
end


function [transmitters_genes, transmitters] = getTransmitorsGenes(path_genes_mat, selectedProbesData)
    path_genes_mat = path_genes_mat';

    transmitters = { 'glutamate' ,'cholinergic','serotonergic','gaba','dopamin'};
    pathwayIds =  [ 4724        ,   4725      ,   4726       , 4727 , 4728    ];
    numberOfTransmitters = length( transmitters );
    numberOfGenes = length( selectedProbesData.entrez_ids );

    transmitters_genes = false(numberOfTransmitters, numberOfGenes);
    for i =1:numberOfTransmitters
        currentPathID = pathwayIds(i);
        transmitors_entrez = find(path_g