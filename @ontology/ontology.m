classdef ontology 
   
    properties
        dependencyMatrix;
        structureLabels;
        structureColors;
        recursiveStructure;
        directedDistanceMatrix;
        unDirectedDistanceMatrix;
    end
    
    methods
        function undirectedMatrix = getUndirectedMatrix(obj)
            undirectedMatrix = obj.dependencyMatrix + obj.dependencyMatrix';
        end
        
        
        function childOfNodes = getIndexesOfChilds(obj,parentNodeIndexes)
            [allChilds, ~] = obj.allChildNodes();
            childOfNodes = allChilds(parentNodeIndexes, :);
            childOfNodes = sum(childOfNodes,1);
        end

        function [allChilds, nodeLevel] = allChildNodes(obj)
            allChilds = inv(eye(size(obj.dependencyMatrix)) - obj.dependencyMatrix);
            nodeLevel = sum(allChilds,1);
        end
        
        
        function reducedOntology = reduceToLeafAndParents(obj,leafIndices)

            [allChilds, ~] = obj.allChildNodes();
            leafParents = allChilds(:,leafIndices);
            leafsAndParents = any(leafParents,2);
            reducedOntology = obj.reduceOntologyToIndexes(leafsAndParents);
        end
        
        function reducedOntology = reduceToNodeAndChilds(obj,nodesIndices)

            [allChilds, ~] = obj.allChildNodes();
            nodeChilds = allChilds(nodesIndices,:);
            nodeAndChilds = any(nodeChilds,1);
            reducedOntology = obj.reduceOntologyToIndexes(nodeAndChilds);
        end
        
        function reducedOntology = reduceOntologyToIndexes(obj,appears)
            reducedOntology = ontology;
            reducedOntology.dependencyMatrix = obj.dependencyMatrix(appears,appears);
            reducedOntology.structureLabels = obj.structureLabels(appears,:);
            reducedOntology.unDirectedDistanceMatrix = obj.unDirectedDistanceMatrix(appears,appears);
            reducedOntology.directedDistanceMatrix = obj.directedDistanceMatrix(appears,appears);
        %    reducedOntology.adjacancyMatrix = ontology.adjacancyMatrix(appears,appears);
        %    reducedOntology.bellowThresholdIndices = ontology.bellowThresholdIndices(appears,appears);
        %    reducedOntology.similarityMeasure = ontology.similarityMeasure(appears,appears);

        end
    end
    
end