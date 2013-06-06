function colors = createColorMap(numOfElements)
    colormapValues = colormap(jet);
    indexOfColor = round(1:size(colormapValues,1)/numOfElements:size(colormapValues,1));
    colors = colormapValues(indexOfColor,:);
end