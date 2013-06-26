function newColorMap = createJetColorMap(numberOfElements, colorMap)
% start at blue 0,0,0.5625
% go to light blue 0,0,1
% then add to the green chanel
% then add to the red and reduce the blue
% then reduce the green
% then reduce the red until 0.5,0,0

    if ~exist('colorMap','var')
        colorMap = jet;
    end

    numElementsInMap = size(colorMap,1);
    defaultTicks = 1:numElementsInMap;

    x = 1 : (numElementsInMap-1)/numberOfElements : numElementsInMap;
    newColorMap = interp1(defaultTicks,colorMap,x,'linear');

end