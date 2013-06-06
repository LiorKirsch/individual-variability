% load('pcaDataAllRegions.mat');
% 
% numberOfSamplesPerArea = sum(experimentsLocationMatrix,1);
% load('humanOntologyObject');


hasAboveKButAllChildsDoesNot = getNodesWithAtLeastKSamples(experimentsLocationMatrix, 10);

subSetIndices = find(hasAboveKButAllChildsDoesNot);
indicesInSubset = subset2ChildTable' * subSetIndices ;
indicesOfExperiments = experimentsLocationMatrix * indicesInSubset;

childOfNodes = humanOntology.getIndexesOfChilds(2);
[allChilds,nodeLevel] = humanOntology.allChildNodes();

 nubmerOfSamplesInNode = zeros(size(allChilds,1),1);
for i=1:size(allChilds,1)
    childsOfNode = logical(allChilds(i,:));
    experimentsOfNode = experimentsLocationMatrix(:,childsOfNode);
    nubmerOfSamplesInNode(i) = sum(sum(experimentsOfNode));
    
end
    

k = 10;
isEndNode = nan(size(allChilds,1),1);
hasAboveKSamples = nubmerOfSamplesInNode >= k;
allChildsHaveMoreThenK = false(size(allChilds,1),1);
allChildsHaveLessThenK = false(size(allChilds,1),1);
for i=1:size(allChilds,1)
    onlyDirectChilds = logical(full(humanOntology.dependencyMatrix(i,:)));
    childNodeSamples = nubmerOfSamplesInNode(onlyDirectChilds);
    allChildsHaveMoreThenK(i) = all(childNodeSamples >= k);
    allChildsHaveLessThenK(i) = all(childNodeSamples < k);
end

hasAboveKButChildsDoesNot = hasAboveKSamples & (~allChildsHaveMoreThenK);
hasAboveKButAllChildsDoesNot = hasAboveKSamples & allChildsHaveLessThenK;

%some checks
onlyChilds = allChilds - eye(size(allChilds));
foundAlsoInChild = false;
for i=1:size(onlyChilds,1)
    if hasAboveKButAllChildsDoesNot(i)
       childsOfNode = logical(onlyChilds(i,:))';
        foundAlsoInChild = any( childsOfNode & hasAboveKButAllChildsDoesNot);
    end
end