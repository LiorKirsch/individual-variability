function normalizedMatrix = normalizeMatrixByDim(matrixToNormalize, dim)

    sizeOfDim = size(matrixToNormalize);
    meanOfMatrix = mean(matrixToNormalize,dim);
    sizeToRepeat = sizeOfDim ./ size(meanOfMatrix);
    
    normalizedMatrix = (matrixToNormalize - repmat( meanOfMatrix, sizeToRepeat ) ) ./ repmat( std(matrixToNormalize,1, dim), sizeToRepeat );

end