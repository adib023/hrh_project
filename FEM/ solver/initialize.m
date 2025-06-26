function [nDof,eqNum, mat_stiff,Wmat, forceNode, resNode] =  ... 
			initialize(nNode,xNode,pbcNode,matProp,leftNode, rightNode, ...
            frcCoord_and_fctr, res_anal_node_coord)

global nMatType DIM 

presNode = [leftNode; rightNode];

nDim = DIM ;
nDof = 0;
eqNum = zeros(nNode,3) ;

for p = 1:nNode
    for q = 1:nDim

        if (eqNum(p,q) == 0)
            nDof = nDof + 1 ;
            eqNum(p,q) = nDof ;
        end
    end
end



nPBC = size(pbcNode,1);

for mat =  1:nMatType

    E = matProp(mat,1);
    nu = matProp(mat,2);
    rho = matProp(mat,3);

    mat_stiff(mat,1) = E/(2*(1+nu));
    mat_stiff(mat,2) = E*nu/(1+nu)/(1-2*nu);
end

% constructing Wmat 
masterNodes = pbcNode(:,1) ; 
followerNodes = pbcNode(:,2) ; 

colID = 0 ; 
count  = 0 ; 

for p = 1:nNode
	for q = 1:nDim 

		rowID = nDim*(p-1) + q ;  

        if (ismember(p,followerNodes) == false) &&...
            (ismember(p,presNode) == false)
            	colID = colID + 1 ;

    			count = count + 1 ;

    			rowVec(count) = rowID ;
    			colVec(count) = colID ;
    			val(count) = 1.0 ;
            
        end

		if (ismember(p,masterNodes) == true)

			indx = followerNodes( find(masterNodes == p) ) ;  
			
			rowID = nDim*(indx-1) + q ; 

			count = count + 1 ; 

			rowVec(count) = rowID ; 
			colVec(count) = colID ; 
			val(count) = 1.0 ; 
		end
	end
end

Wmat = sparse( rowVec , colVec , val , nNode*nDim , colID ) ;  


%%%% find force application nodes 

for p = 1:size(frcCoord_and_fctr, 1)
    targetCoord = frcCoord_and_fctr(p,1:3);
    checkNorm = sqrt(sum((xNode -  targetCoord).^2,2));
    % FroceNodeID =  find(checkNorm < 10^-6);
    %checkMin = find(checkNorm == min(checkNorm))
    % xNode(checkMin,:)
    % p
    % FroceNodeID
    FroceNodeID =  find(checkNorm == min(checkNorm));

    forceNode(p,1) = FroceNodeID;
    forceNode(p,2:4) =  frcCoord_and_fctr(p,4:6);
end

%%% determine the indx of res_coord_list


for p = 1:size(res_anal_node_coord, 1)
    targetCoord = res_anal_node_coord(p,1:3);
    checkNorm = sqrt(sum((xNode -  targetCoord).^2,2));
    resNodeID =  find(checkNorm == min(checkNorm));
    resNode(p) = resNodeID;
end


end
