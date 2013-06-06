function analyzeGoCat()
    %load('catScores.mat');
    load('catScoresNotRemove.mat');
    

    meanCatScores = squeeze(mean(catDistanceScore,2));
    
    plotScoresAndCatSize(numOfGenesInCategory, meanCatScores,regionNames);
  
    [a,b,c] = unique(numOfGenesInCategory);
    %subsampleAndCheckPrecision(a, 'newRandomClasses.mat');
   
   
    random1 = load('randomSizePrecision.mat','percentInCategory','sizesOfClasses');
%     random2 = load('newRandomClasses.mat','percentInCategory','sizesOfClasses');
    
    
    
%     randomScores = cat(2, random1.percentInCategory, random2.percentInCategory);
%     catSizes = cat(1, random1.sizesOfClasses, random2.sizesOfClasses);
    randomScores = random1.percentInCategory;
    catSizes = random1.sizesOfClasses;
    
    meanRandomScores = mean(randomScores,2);
    meanCatScores = mean(catDistanceScore,2);
    
    
    pvalues = nan(length(cat_ids),1);
    decisions = nan(length(cat_ids),1);
    
    [numOfGenesInCategory ,sortIndecies]= sort(numOfGenesInCategory);
    catDistanceScore = catDistanceScore(sortIndecies);
    cat_ids = cat_ids(sortIndecies);
    go_gene_mat = go_gene_mat(sortIndecies,:);
    meanCatScores = meanCatScores(sortIndecies);
    
    for i = 1:length(cat_ids)
        numberOfGenesInCurrentCat = numOfGenesInCategory(i);
        randomSampleIndex = catSizes == numberOfGenesInCurrentCat;
        assert( (sum(randomSampleIndex)==1) );
        randomScores = meanRandomScores(:,randomSampleIndex);
        [pvalues(i),decisions(i),stats] = ranksum(randomScores,meanCatScores(i) ) ;
    end
    
end

function plotScoresAndCatSize(catSize, catScores, regionNames)

    minimumValue =  min(min(catScores));
    [uniqueValues,indices, reverseIndices] = unique(catSize);
    
    for i = 1:length(regionNames)
        
       
        h  = figure(i*23);
        hold on;
        %subplot(1,2,1);
        plot(catSize,catScores(:,i),'.');
        ylabel('precision at 1','fontSize',20); 
        xlabel('category size','fontSize',20);
        title(regionNames{i},'fontSize',20);
        %ylim([(minimumValue -0.05) 1.05]);
        
        
        %subplot(1,2,2);
        
%         meanValues = nan(size(uniqueValues));
%         stdValues = nan(size(uniqueValues));
%         for j = 1:length(uniqueValues)
%             currentSizeIndices = reverseIndices == j;
%             meanValues(j) = mean(catScores(currentSizeIndices,i));
%             stdValues(j) = std(catScores(currentSizeIndices,i));
%         end
%         errorbar(uniqueValues, meanValues, stdValues,'.r');
%         

        set(gca,'xscale','log')
        fileName = ['regionsFigures\GoPrecision',regionNames{i},'.png'];
        saveas(h,fileName) ;
        hold off;
    end

end