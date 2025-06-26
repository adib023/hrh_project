function writePolyBasis(nNode,xNode,outputFile)

% routine will write polynomial basis to a file: use for debugging 
% and performance comparison  

% steps: 
% 1. find x_max , y_max, z_max . 
% 2. find mode values at each node using the following local coords: 
%	alpha = x/x_max ; beta = y/y_max ; gama = z/z_max  
% 3. write mode to file
%

nMode = 30 ; 

% step 1 : 
xmin 	= 	min(xNode(:,1)) ; 
xmax 	= 	max(xNode(:,1)) ; 
Lx 	=	xmax - xmin ;

ymin 	=	min(xNode(:,2)) ; 
ymax 	= 	max(xNode(:,2)) ; 
Ly 	=	ymax - ymin ;

zmin 	= 	min(xNode(:,3)) ; 
zmax 	= 	max(xNode(:,3)) ; 
Lz 	=	zmax - zmin ;

% step 2:
for p = 1:nNode

	alph = ( xNode(p,1) - xmin ) / Lx ;  
	beta = ( xNode(p,2) - ymin ) / Ly ;  
	gama = ( xNode(p,3) - zmin ) / Lz ;  

	modeX(p,1:3,1:nMode) = getPolyMode(alph,beta,gama) ; 
end


% step 3: 
fid = fopen(outputFile,'w') ;

fprintf(fid,'%d\n',nMode);
for p = 1:nMode
	for q = 1:3 % displacement 
		fprintf(fid,'%.16f\t',modeX(:,q,p)) ; 
		fprintf(fid,'\n') ; 
	end
	for q = 1:3 % velocity
		fprintf(fid,'%.16f\t',modeX(:,q,p)) ; 
		fprintf(fid,'\n') ; 
	end
		fprintf(fid,'\n') ; 
end

end

function [modeVal] = getPolyMode(alph,beta,gama)

modeVal = zeros( 3 , 30 ) ; 

mvec(1) = 1.0 ; 

mvec(2) = alph ; 
mvec(3) = beta ; 
mvec(4) = gama ; 

mvec(5) = alph * alph ; 
mvec(6) = beta * beta ; 
mvec(7) = gama * gama ; 

mvec(8)  = beta * gama ; 
mvec(9)  = gama * alph ; 
mvec(10) = alph * beta ; 

modeVal(1,1:10) = mvec(1:10) ; 
modeVal(2,11:20) = mvec(1:10) ; 
modeVal(3,21:30) = mvec(1:10) ; 

end

