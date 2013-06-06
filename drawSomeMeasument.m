function drawSomeMeasument(experimentsDataMatrix, experimentsSubjectMatrix, smoothParam)
    if ~exist('smoothParam')
        smoothParam = 1;
    end
   
    figure(2344212);
    meanForGene = nan(6, size(experimentsDataMatrix,2));
    dataPersoni = experimentsDataMatrix(experimentsSubjectMatrix == 1,:);
    [meanForGene(1,:), sortOrder] = sort(mean(dataPersoni,1));
    
    for i = 1:6
        dataPersoni = experimentsDataMatrix(experimentsSubjectMatrix == i,:);
        meanForGeneTemp = mean(dataPersoni,1);
        meanForGene(i,:) = smooth(meanForGeneTemp(sortOrder),smoothParam);
    end
    meanMean = std(meanForGene,0,2);
    [sortedMean, sortedByMeanMean] = sort(meanMean,'descend');
    plot(meanForGene(sortedByMeanMean,:)','.');
    
    xlabel('gene index (sorted by person 1 scores)');
    ylabel('mean (across samples)');
    
    legend([repmat('human',[6,1]), num2str(sortedByMeanMean)]); 
    
     figure(234422);
    stdForGene = nan(6, size(experimentsDataMatrix,2));
    dataPersoni = experimentsDataMatrix(experimentsSubjectMatrix == 1,:);
    %[stdForGene(4,:), sortOrder] = sort(std(dataPersoni,1));
    
    for i =1:6
        dataPersoni = experimentsDataMatrix(experimentsSubjectMatrix == i,:);
        stdForGeneTemp = std(dataPersoni,0,1);
        stdForGene(i,:) = smooth(stdForGeneTemp(sortOrder),smoothParam);
    end
    stdStd = std(stdForGene,0,2);
    [sortedStd, sortedByMeanStd] = sort(stdStd,'descend');
    plot(stdForGene(sortedByMeanStd,:)','.');
    
    xlabel('gene index (sorted by person 1 scores)');
    ylabel('std (across samples)');
    
    legend([repmat('human',[6,1]), num2str(sortedByMeanStd)]); 
    
    
    
end

