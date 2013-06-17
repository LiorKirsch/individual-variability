function analyzeGoCat()
    load('catScores.mat');
    %load('catScoresNotRemove.mat');
    
    load('randomSizePrecision','percentInCategoryUsingGenes', 'percentInCategory','sizesOfClasses');
    houseKeeping = load('catHouseKeepingScores.mat', 'catDistanceScore', 'onlyCatDistanceScore','cat_ids','numOfGenesInCategory');
    
    subsetConf.onlyBrain = false;
    subsetConf.minCatSize = 50;
    subsetConf.maxCatSize = 300 ;
    subsetConf.catBranch = 'P';
    withOutGenes = true;
    
    [cat_ids, aspects, go_gene_mat, numOfGenesInCategory, catDistanceScore, onlyCatDistanceScore] = chooseSubsetOfCategories(subsetConf, cat_ids, aspects, go_gene_mat, numOfGenesInCategory, catDistanceScore, onlyCatDistanceScore);

    cat_ids = [1, cat_ids];
    numOfGenesInCategory = [560 ; numOfGenesInCategory];
    catDistanceScore = cat(1, houseKeeping.catDistanceScore, catDistanceScore);
    onlyCatDistanceScore = cat(1, houseKeeping.onlyCatDistanceScore, onlyCatDistanceScore);
    
    if withOutGenes
        meanRandCatScores = squeeze(mean(percentInCategory,2));
        meanCatScores = squeeze(mean(catDistanceScore,2));
    else
        catDistanceScore = onlyCatDistanceScore;
        meanRandCatScores = squeeze(mean(percentInCategoryUsingGenes,2));
        meanCatScores = squeeze(mean(onlyCatDistanceScore,2));
    end
    
    
    
    
   [sigScoresAbove, sigScoresBelow] = getSignificeNormalScores(numOfGenesInCategory, meanCatScores,meanRandCatScores,full(sizesOfClasses),regionNames);
   correctedNormalScoresAbove = decideSignificance(sigScoresAbove);
   correctedNormalScoresBelow = decideSignificance(sigScoresBelow);
%    sigScores = getSignificeScores(numOfGenesInCategory, meanCatScores,meanRandCatScores,full(sizesOfClasses),regionNames);
    

    subSetNames = {'Cerebellar Cortex', 'Cingulate gyrus', 'Occipital Lobe'};
    %subSetNames =regionNames;
    
    for i =1:length(subSetNames)
       subsestResultsAbove{i} =  correctedNormalScoresAbove(:, strcmp(regionNames,subSetNames{i}) );
       subsestResultsBelow{i} =  correctedNormalScoresBelow(:, strcmp(regionNames,subSetNames{i}) );
    end
   
    if subsetConf.onlyBrain
        analyzeCatWithScores(cat_ids, aspects, subsestResultsAbove, subSetNames, 'variBrainAboveMean.html','without cat variability among subjects is stronger',numOfGenesInCategory);
        analyzeCatWithScores(cat_ids, aspects, subsestResultsBelow, subSetNames, 'variBrainBelowMean.html','without cat variability is reduced', numOfGenesInCategory);
    else
        analyzeCatWithScores(cat_ids, aspects, subsestResultsAbove, subSetNames, 'variAboveMean.html','without cat variability among subjects is stronger',numOfGenesInCategory);
        analyzeCatWithScores(cat_ids, aspects, subsestResultsBelow, subSetNames, 'variBelowMean.html','without cat variability is reduced', numOfGenesInCategory);
    end
   
    plotScoresAndCatSize(numOfGenesInCategory, meanCatScores,meanRandCatScores,sizesOfClasses,regionNames);
    
%     meanOnlyCatScores = squeeze(mean(onlyCatDistanceScore,2));
%     plotScoresAndCatSize(numOfGenesInCategory, meanOnlyCatScores,regionNames);
    
  
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

function correctedScores = decideSignificance(scores)
    numberOfGenes = size(scores,1);
    numberOfRegions = size(scores,2);
        
    newScores = reshape (scores, [numel(scores),1]);
    correctedScores = mafdr(newScores, 'BHFDR', true);
    correctedScores = reshape(correctedScores, size(scores));
   
end
function sigScores = getSignificeScores(catSize, catScores,randCatScores,sizesOfClasses,regionNames)
    sigScores = nan(size(catScores));
    for regionIter =1:length(regionNames)
        for catIter =1:size(catScores,1)
            sizeIndex = sizesOfClasses == catSize(catIter);
            randomSamples = squeeze(randCatScores(:,regionIter,sizeIndex));
            sigScores(catIter,regionIter)  = ranksum(randomSamples , catScores(catIter,:));
        end
    end
end


function [cat_ids, aspects, go_gene_mat, numOfGenesInCategory, catDistanceScore, onlyCatDistanceScore] = chooseSubsetOfCategories(subsetConf, cat_ids, aspects, go_gene_mat, numOfGenesInCategory, catDistanceScore, onlyCatDistanceScore)
    
    subsetIndices = true(length(cat_ids),1);
    
    % only Biological processes cat
      biologicalProcessesCat  = strcmp(aspects,subsetConf.catBranch);
      subsetIndices = subsetIndices & biologicalProcessesCat';
      
    if subsetConf.onlyBrain
        brainRelated = load('/home/lab/noalis/work3/for_uri_feb_2012/new/old_data/images_go_genes_mat_brain_screened', 'go_cat_names', 'cat_ids');
        onlyBrainRelated = ismember(cat_ids, brainRelated.cat_ids);
        subsetIndices = subsetIndices & onlyBrainRelated';
    end

    midSizeCategories = (numOfGenesInCategory >= subsetConf.minCatSize) & (numOfGenesInCategory <= subsetConf.maxCatSize);
    subsetIndices = subsetIndices & midSizeCategories;
    
    aspects = aspects(subsetIndices);
    cat_ids = cat_ids(subsetIndices);
    catDistanceScore = catDistanceScore(subsetIndices,:,:);
    onlyCatDistanceScore = onlyCatDistanceScore(subsetIndices,:,:);
    go_gene_mat = go_gene_mat(subsetIndices ,:);
    numOfGenesInCategory = numOfGenesInCategory(subsetIndices);
    
end
function [sigScoresAbove, sigScoresBelow] = getSignificeNormalScores(catSize, catScores,randCatScores,sizesOfClasses,regionNames)
    sigScoresAbove = nan(size(catScores));
    sigScoresBelow = nan(size(catScores));
    for regionIter =1:length(regionNames)
        for catIter =1:size(catScores,1)
            sizeIndex = sizesOfClasses == catSize(catIter);
            randomSamples = squeeze(randCatScores(:,regionIter,sizeIndex));
            [~, sigScoresAbove(catIter,regionIter)] = ztest(catScores(catIter,regionIter),mean(randomSamples),std(randomSamples), 0.05, 'right');
            [~, sigScoresBelow(catIter,regionIter)] = ztest(catScores(catIter,regionIter),mean(randomSamples),std(randomSamples), 0.05, 'left');
        end
    end
end

function plotScoresAndCatSize(catSize, catScores,meanRandCatScores, sizeOfClasses,regionNames)

    minimumValue =  min(min(catScores));
    maximumValue =  max(max(catScores));
    [uniqueValues,indices, reverseIndices] = unique(catSize);
    
    meanRandomScores = squeeze(mean(meanRandCatScores,1));
    stdRandomScores = squeeze(std(meanRandCatScores,1,1));
    for i = 1:length(regionNames)
        
       
        h  = figure(i*23);
        hold on;
        %subplot(1,2,1);
        plot(catSize,catScores(:,i),'.');
        
        errorbar(sizeOfClasses, meanRandomScores(i,:), stdRandomScores(i,:), 'xr');
        
        ylabel('between / within','fontSize',20); 
        xlabel('category size','fontSize',20);
        title(regionNames{i},'fontSize',20);
        %ylim([(minimumValue -0.1) (maximumValue +0.1)]);
        
        
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