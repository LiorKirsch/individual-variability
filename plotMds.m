function mdsEmbeding = plotMds(experimentsDataMatrix,experimentRegion, experimentsSubjectMatrix, regionNames,regionColors)

    expressionDistances = distanceUsingMatrices(experimentsDataMatrix);
%     all(all(expressionDistances == expressionDistances'))
%     all(all(expressionDistances >= 0 ))
%     diagIndices = 1:size(expressionDistances,1)+1:numel(expressionDistances); 
%     all(expressionDistances(diagIndices) == 0 )
    
    pdistExpressionDistance = pdist(experimentsDataMatrix,'euclidean');
    expressionDistances2 = squareform(pdistExpressionDistance);
    expressionDistances2 = expressionDistances2.^2;
    all(expressionDistances2(:) == expressionDistances(:))
    
%     [cmdsEmbeding,eigvals] = cmdscale(pdistExpressionDistance);
    [mdsEmbeding,stress] = mdscale(pdistExpressionDistance,2,'criterion','metricstress');
%     Y1 = mdscale(pdistExpressionDistance,2);
%     Y2 = mdscale(pdistExpressionDistance,3);


    numberOfRegions = size(experimentRegion,2);
    experimentRegionIndex = experimentRegion *  (1:numberOfRegions)';
    experimentRegionColor = regionColors(experimentRegionIndex,:);
    
    figure(1);
    scatterRegionColor(mdsEmbeding, experimentRegion, experimentRegionColor,regionNames);
    
    
    numberOfPeople = 6;
    experimentsSubjectMatrixLogical = experimentsSubjectMatrix;
    for i =1:numberOfPeople
        subjectNames{i} = sprintf('human%d',i);
    end
    
    subjectColor = createColorMap(6);
    experimentSubjectColor = experimentsSubjectMatrixLogical * subjectColor;
    figure(2);
    scatterRegionColor(mdsEmbeding, experimentsSubjectMatrixLogical, experimentSubjectColor,subjectNames);
    
    for i =1:numberOfRegions
        h = figure(i*321);
        onlyInregion = experimentRegion(:,i);
        scatterRegionColor(mdsEmbeding(onlyInregion,:), experimentsSubjectMatrixLogical(onlyInregion,:), experimentSubjectColor(onlyInregion,:),subjectNames);
        title(regionNames(i));
        ylim([-300 300]);
        xlim([-300 300]);
        fileName = ['regionsFigures\mds-',regionNames{i},'.png'];
        saveas(h,fileName) ;
    end
end


function scatterRegionColor(score, experimentCategory, experimentColors,categoryNames)
    numberOfCategories = size(experimentCategory,2);
%     experimentCategoryIndex = experimentCategory *  (1:numberOfCategories)';
% %     experimentCategoryLabel = categoryNames(experimentCategoryIndex,:);
%     experimentCategoryColor = categoryColors(experimentCategoryIndex,:);
    markerType = {'o' ,'^','v','x','s','d', '>','+','p' ,'h','s' ,'<','x','o','d'};

    
    for j=1: numberOfCategories
       tmpIndex = experimentCategory(:,j);
%        tmpRegion =  experimentCategoryLabel(tmpIndex);
       tmpColor = experimentColors(tmpIndex,:);
       tmpIExperiments = score(tmpIndex,:);
       
       hold on;
%        assert( all(tmpColor(1,:) == experimentColors(j,:) ) );
       scatter(tmpIExperiments(:,1),tmpIExperiments(:,2),10,tmpColor(1,:),markerType{j});

%         figure(2);
%         hold on;
%         scatter(tmpIExperiments(:,2),tmpIExperiments(:,3),10,tmpColor(1,:),markerType{j});
%         figure(3);
%         hold on;
%         scatter(tmpIExperiments(:,3),tmpIExperiments(:,4),10,tmpColor(1,:),markerType{j});

    end
    xlabel('MDS #1');     ylabel('MDS #2');
    legend(categoryNames);
    figureHandle = gcf;     set(findall(figureHandle,'type','text'),'fontSize',14);
%     figure(2);
%     xlabel('MDS #2');     ylabel('MDS #3');
%     legend(categoryNames);
%     figureHandle = gcf;     set(findall(figureHandle,'type','text'),'fontSize',14);
%     figure(3);
%     xlabel('MDS #3');     ylabel('MDS #4');
%     legend(categoryNames);
%     figureHandle = gcf;     set(findall(figureHandle,'type','text'),'fontSize',14);
end

function squaredDistances = distanceUsingMatrices(x)
    Qx=repmat(dot(x,x,2),1,size(x,1));
    squaredDistances = (Qx+Qx'-2*(x*x'));
end