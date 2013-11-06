function writeCellToFile(fileName, cellArray)

    fprintf('writing to %s\n', fileName);
    fid = fopen(fileName,'w');
    
    for i  = 1:length(cellArray)
       fprintf(fid, '%s\n', cellArray{i});
    end

end