function prepareAgeData()

    [age_pEnsmble, age_gEnsmble, age_ageNumber, age_ageDescription] = textread('geneAge.csv', '%s %s %d %s', 'delimiter' , ',', 'headerlines' ,2);
    [mart_EnsemblGeneID, mart_EntrezGeneID,HGNC_symbol] = textread('martExport2.csv', '%s %d %s', 'delimiter' , ',', 'headerlines' ,1);

    easyData = load('easyFormatHumanData.mat');

    selectedProbesData = easyData.selectedProbesData;

    ages = zeros(length(selectedProbesData.entrez_ids),19);
    [~, someIndexes] = ismember( 1:19, age_ageNumber);
    agesDescription = age_ageDescription( someIndexes);

    for i = 1:length(selectedProbesData.entrez_ids)
       geneEntrez = selectedProbesData.entrez_ids(i);
       entrezRecord = find (mart_EntrezGeneID == geneEntrez);

       for j = 1:length(entrezRecord)
          ensmbleId =  mart_EnsemblGeneID( entrezRecord(j) );
          ageIndex = strcmp(ensmbleId,age_gEnsmble);
          ages(i, age_ageNumber(ageIndex)) = 1;
       end
    end

    hist(sum(ages,2));
    easyData.ages = ages;
    easyData.agesDescription = agesDescription;


    save('easyFormatHumanData_withAges.mat','-struct', 'easyData');
end