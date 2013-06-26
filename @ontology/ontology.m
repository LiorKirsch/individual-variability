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
        
        function ontologyColors = createColors(obj)
           ontologySize = size(obj.dependencyMatrix,1);
           newColorMap = createJetColorMap(ontologySize);
           colors = ones(ontologySize,1);

%            colors = obj.recursiveColor(2, 0,1,colors);
           
           segments = linspace(0,1,9) ;
           colors = obj.recursiveColor(4, segments(1),segments(end-5),colors);
           colors = obj.recursiveColor(302, segments(end-5),segments(end-4),colors);
           colors = obj.recursiveColor(423, segments(end-4),segments(end-3),colors);
           colors = obj.recursiveColor(752, segments(end-3),segments(end-2),colors);
           colors = obj.recursiveColor(923, segments(end-2),segments(end-1),colors);
           colors = obj.recursiveColor(1114, segments(end-1),segments(end),colors);
           
           colorIndex = ceil( colors * ontologySize);
            ontologyColors = newColorMap(colorIndex,:);
        end
        function colors = recursiveColor(obj, nodeIndex, startInterval,endInterval,colors)
            meanValue = mean([startInterval,endInterval]);
            colors(nodeIndex) = meanValue;
            childOfNodes = find(obj.dependencyMatrix(nodeIndex,:));
            
            segments = linspace(startInterval,endInterval,length(childOfNodes)+1) ;
            for i=1: length(childOfNodes)
               childNodeIndex = childOfNodes(i);
               colors =  obj.recursiveColor(childNodeIndex, segments(i),segments(i+1),colors);
            end
        end
        function regionColors = getColorByRegionName(obj, regionNames)
            colors = obj.createColors();
            [~,indexInOntlogy] = ismember(regionNames, obj.structureLabels(:,4));
            regionColors = colors(indexInOntlogy,:); 
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