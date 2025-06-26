function [dV1g, dV2g, dV3g, dV4g]= magnet_dvxcmdu_assmebly(shapeFunction, nGaussPts, gaussPts, gaussWeights,...
    indx,xNode,UUR,elemMagnetID)

global DIM;
numMagnet = 4;
numNodeElem = size(indx,2);

func =@plus;
% dVg_cell = cell(4,1);

% dVg_cell{1} =  sparse(DIM,DIM*size(xNode,1));
% dVg_cell{2} =  sparse(DIM,DIM*size(xNode,1));
% dVg_cell{3} =  sparse(DIM,DIM*size(xNode,1));
% dVg_cell{4} =  sparse(DIM,DIM*size(xNode,1));

dVg_whole  = sparse(numMagnet*DIM,DIM*size(xNode,1));

for elem = 1:size(indx,1)
    matID = elemMagnetID(elem);
    if matID == 0
        continue
    end
    elemNode = indx(elem,:);

    Xel = xNode(elemNode,:)';
    Uel = UUR(elemNode,:)';

	[Xcm] = getElemCM(DIM,numNodeElem,Xel,shapeFunction, ...
			nGaussPts , gaussPts , gaussWeights) ;

    dVxcmdU =  zeros(DIM,DIM*numNodeElem);

    parfor p = 1:nGaussPts
        r =  gaussPts(p,1:DIM);
        w = gaussWeights(p);

        [N,DN] = shapeFunction(r);
        J = Xel*DN;
        detJ =  det(J);

        invJ = inv(J);
        B = (DN*invJ);
        NoI = kron(N',eye(3))';
        defF = Uel*B + eye(3);
        deterF =  det(defF);
        BoI = kron(B',eye(3));
        cofF = (deterF*inv(defF))';
        cofF = cofF(:); %%%%% 9 by 1

        dVxcmdU = dVxcmdU + (NoI*deterF + (Xcm + Uel*N')*(BoI'*cofF)')*detJ*w;
    end
    
    % dV_elem_cell{matID} = dVxcmdU; 
    % dVg_whole = assemble_dVxcmdU(DIM,size(xNode,1),elemNode,matID, dVxcmdU, dummy_rowIndxMat); %%%%% why this error is here?????

    [rowIndxMat, colIndxMat] =  get_colIndxMat(DIM, elemNode, matID, numNodeElem);
    dVg_whole= func(dVg_whole, sparse(rowIndxMat,colIndxMat,dVxcmdU,numMagnet*DIM,DIM*size(xNode,1)));
    
end

dV1g = dVg_whole(1:3,:);
dV2g = dVg_whole(4:6,:);
dV3g = dVg_whole(7:9,:);
dV4g = dVg_whole(10:12,:);


end

% function dVg_whole = assemble_dVxcmdU(DIM, nNode, elemNode, matID, dVg_whole, dVxcmdU, dummy_rowIndxMat)
% 
% fun = @plus;
% colIndx = [];
% 
% for n = 1:length(elemNode)
%     node = elemNode(n);
% 
% 	lim = (DIM*node - DIM + 1) : ( DIM*node ) ; 
% 	colIndx = [colIndx lim];   
% end
% 
% rowIndxMat =  dummy_rowIndxMat + DIM*(matID-1);
% colIndxMat =  [colIndx;colIndx;colIndx];
% 
% dVg_cell{matID} = fun(dVg_cell{matID},sparse(rowIndxMat,colIndxMat,dVxcmdU,DIM, DIM*nNode)); 
% 
% end

function [rowIndxMat,colIndxMat] = get_colIndxMat(DIM, elemNode, matID, numNodeElem)

colIndx = [];

for n = 1:length(elemNode)
    node = elemNode(n);

	lim = (DIM*node - DIM + 1) : ( DIM*node ) ; 
	colIndx = [colIndx lim];   
end

rowIndxMat =  [ones(1,3*numNodeElem); ones(1,3*numNodeElem)*2; ones(1,3*numNodeElem)*3] + DIM*(matID - 1);
colIndxMat =  [colIndx;colIndx;colIndx];


end





% function dVg_cell = assemble_dVxcmdU(DIM, elemNode, matID, dVg_cell, dVxcmdU )
% % dV1 = dV_cell{1};
% % dV2 = dV_cell{2};
% % dV3 = dV_cell{3};
% % dV4 = dV_cell{4};
% 
% for n = 1:length(elemNode)
%     node = elemNode(n);
% 
% 	lim = (DIM*node - DIM + 1) : ( DIM*node ) ; 
% 	loc_lim = (DIM*n-DIM+1) : (DIM*n) ; 
% 
%     % dV1g(:,lim) = dV1g(:,lim) + dV1(:,loc_lim);
%     % dV2g(:,lim) = dV2g(:,lim) + dV2(:,loc_lim);
%     % dV3g(:,lim) = dV3g(:,lim) + dV3(:,loc_lim);
%     % dV4g(:,lim) = dV4g(:,lim) + dV4(:,loc_lim);
% 
%     dummy_dVxcmdU =  zeros(DIM,DIM*size(xNode,1));
%     dummy_dVxcmdU(:,lim) =  dVxcmdU(:,loc_lim);
% 
% end
% 
% dV_cell{matID} = dV_cell{matID} + dummy_dVxcmdU; 
% 
% end



