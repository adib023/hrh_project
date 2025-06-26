function [Xcm,vol] = getElemCM(DIM,Xel, shapeFunction, ...
			nGaussPts , gaussPts , gaussWeights)

% routine to get center of mass of element 
vol = 0 ; %xvol = zeros(DIM,numNodeElem) ; 
xvol = zeros(DIM,1);
for p = 1:nGaussPts 

	r = gaussPts(p,1:DIM) ; 
	w = gaussWeights(p) ; 

	[N,DN] = shapeFunction(r);
    Nx =  Xel * N';

        J = Xel*DN;
       	detJ =  det(J);

	xvol = xvol + Nx*detJ*w ;  
	vol  = vol + detJ*w ; 
end

Xcm	=	xvol / vol ; 


end
