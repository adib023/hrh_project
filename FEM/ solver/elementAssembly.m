function [rowIndxVec,colIndxVec, fIndxVec] = elementAssembly(probType,elemNode, ...
			eq_num )

nNodeElem 	= 	size(elemNode,2) ;
nDim		=	3 ;
Mvec = [];

count = 0 ;
fcount = 0 ; 
for p = 1:nNodeElem

	node_p = elemNode(p) ;
	for px = 1:nDim
        % node_p
        % size(eq_num)
		rowID = eq_num(node_p,px) ; 

%		Pfull(rowID) = Pfull(rowID) + Pelem(nDim*(p-1)+px) ; 

	fcount = fcount + 1 ; 
	fIndxVec(fcount) = rowID ; 

		for q = 1:nNodeElem
		
			node_q = elemNode(q) ; 
			for qx = 1:nDim 

	colID = eq_num(node_q,qx);

	count = count + 1; 
	rowIndxVec(nDim*(p-1)+px,nDim*(q-1)+qx) = rowID ; 
	colIndxVec(nDim*(p-1)+px,nDim*(q-1)+qx) = colID ; 

%	Kvec(count) 	   = Kelem(nDim*(p-1)+px,nDim*(q-1)+qx ) ; 
%	Mvec(count) 	   = Melem(nDim*(p-1)+px,nDim*(q-1)+qx ) ; 

			end
		end
	end
end


end
