function writeCSV( fileName, stringCellArrays, numericCellArrays,headers)

    fid = fopen(fileName,'w');
    num_rows = length(stringCellArrays{1});
    
    headerString = '';
    for j =  1:length(headers)
           headerString = sprintf('%s,"%s"' , headerString,headers{j});
    end
    headerString = headerString(2:end);

    fprintf(fid,'%s\n', headerString);
    
    for i = 1:num_rows
       numString = '';
       for j = 1:length(numericCellArrays)
           numString = sprintf('%s,%g' , numString,numericCellArrays{j}(i));
       end
       numString = numString(2:end);
       
       
       cellString = '';
       for j = 1:length(stringCellArrays)
           cellString = sprintf('%s,"%s"' , cellString,stringCellArrays{j}{i});
       end
       cellString = cellString(2:end);
       
       fprintf(fid,'%s,%s\n', cellString, numString);
    end
    
    fclose(fid);
end