function someViewingForAnova(anovaScores, releventExprimentsData , groupData)

    notNaN = ~isnan(anovaScores);
    anovaScores = anovaScores(notNaN);
    releventExprimentsData = releventExprimentsData(notNaN,:);
    
    [minValue, minIndex] = min(anovaScores);
    slightMove = groupData + (rand(size(groupData)) - 0.5) /10;

    scatter(slightMove, releventExprimentsData(minIndex,:),'.' );
    
end