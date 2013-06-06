function [hasAboveKButAllChildsDoesNot ,subset2ChildTable]= getNodesWithAtLeastKSamples(experimentsLocationMatrix, k)
    load('humanOntologyObject');

    [allChilds,~] = humanOntology.allChildNodes();

     nubmerOfSamplesInNode = zeros(size(allChilds,1),1);
    for i=1:size(allChilds,1)
        childsOfNode = logical(allChilds(i,:));
        experimentsOfNode = experimentsLocationMatrix(:,childsOfNode);
        nubmerOfSamplesInNode(i) = sum(sum(experimentsOfNode));

    end
    
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

    subset2ChildTable = allChilds(hasAboveKButAllChildsDoesNot,:);
end