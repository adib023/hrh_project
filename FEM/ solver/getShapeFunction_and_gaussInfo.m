function [shapeFunction, nGaussPts, gaussPts, gaussWeights] =  getShapeFunction_and_gaussInfo(elemType)

global DIM;

if strcmp(elemType,'HEX2') == 1

    shapeFunction	= @shapeFnBrick2 ;

elseif strcmp(elemType,'HEX1') == 1

    shapeFunction	= @shapeFnBrick ;

elseif strcmp(elemType,'TET1') == 1

    shapeFunction	= @shapeFnTetra ;

elseif strcmp(elemType,'TET2') == 1

    shapeFunction	= @shapeFnTetra2 ;

else
    elemType
    error('not yet implemented');
end

[nGaussPts,gaussPts,gaussWeights] = GaussQuad(DIM,elemType) ;


end