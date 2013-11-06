classdef ontology 
   
    properties
        dependencyMatrix;
        structureLabels;
        structureColors;
        recursiveStructure;
        directedDistanceMatrix;
        unDirectedDistanceMatrix;
        closestCommonParentIndex;
        meanDistanceToParent;
        shortDistanceToParent;
        longDistanceToParent;
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
        
        function directedDistences = getDistances(obj)
            X = obj.dependencyMatrix;
            directedDistences = ( inv(eye(size(X)) - X) )^2*X;
            directedDistences(directedDistences ==0) = inf;
            directedDistences(logical(eye(size(X)))) = 0;
        end
        
        
        function obj = distanceToCommonParent2(obj)
            [allChilds, nodeLevel] = obj.allChildNodes();
            nodeLevel = nodeLevel';
            distanceToParent = zeros(size(obj.dependencyMatrix));
            obj.closestCommonParentIndex = zeros(size(obj.dependencyMatrix));
            numberOfNodes = size(obj.dependencyMatrix,1);
            
            for i = 1: (numberOfNodes -1)
                ancestorOf_i_j = repmat(allChilds(:,i),1, numberOfNodes ) & allChilds;
                
                nodeReplicated = repmat(nodeLevel,1, numberOfNodes ) ;
                parentLevel = nodeReplicated .* ancestorOf_i_j;
                        
                [~, bestParentIndex] = max(parentLevel, [], 1 );
                distanceToParent(:,i) = obj.unDirectedDistanceMatrix(i,bestParentIndex);
                obj.closestCommonParentIndex(:,i) = bestParentIndex;
            end
            
            obj.longDistanceToParent =  max(distanceToParent, distanceToParent');
            obj.shortDistanceToParent =  min(distanceToParent, distanceToParent');
            obj.meanDistanceToParent = (distanceToParent + distanceToParent')/2;
        end
        
        function [ closestCommonParentIndex, meanDistanceToParent, shortDistanceToParent, longDistanceToParent] = distanceToCommonParent(obj)
            [allChilds, nodeLevel] = obj.allChildNodes();

            
            closestCommonParentIndex = zeros(size(obj.dependencyMatrix));
            meanDistanceToParent = zeros(size(obj.dependencyMatrix));
            shortDistanceToParent = zeros(size(obj.dependencyMatrix));
            longDistanceToParent = zeros(size(obj.dependencyMatrix));

            numberOfNodes = size(obj.dependencyMatrix,1);
            
            
            for i = 1: numberOfNodes
                for j= i+1 : numberOfNodes
                    ancestorOf_i_j = [ allChilds(:,i) , allChilds(:,j)];
                    ancestorOf_i_j = all( ancestorOf_i_j,2);
                    ancestorsIndices = find(ancestorOf_i_j);
                    
                    [~, bestParentmInIndex] = max(nodeLevel (ancestorOf_i_j) );
                    bestParentmInIndex = ancestorsIndices(bestParentmInIndex);
                    
                    parentDistanceFrom_i = obj.directedDistanceMatrix(i,bestParentmInIndex);
                    parentDistanceFrom_j = obj.directedDistanceMatrix(j,bestParentmInIndex);
                    distances = [parentDistanceFrom_i, parentDistanceFrom_j];
                    sortedDist = sort(distances);
                    shortDistanceToParent(i,j) = sortedDist(1);
                    longDistanceToParent(i,j) = sortedDist(2);
                    meanDistanceToParent(i,j) = mean(distances);
                    closestCommonParentIndex(i,j) = bestParentmInIndex;
                    
                end
            end
            
            closestCommonParentIndex = closestCommonParentIndex + closestCommonParentIndex';
            meanDistanceToParent = meanDistanceToParent + meanDistanceToParent';
            shortDistanceToParent = shortDistanceToParent + shortDistanceToParent';
            longDistanceToParent = longDistanceToParent + longDistanceToParent';
            
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