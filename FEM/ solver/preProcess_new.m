function [elemType,nNode, nElem, xNode, indx ,matProp,elemMatIndx,...
    pbcNode, leftNode, rightNode]...
    = preProcess_new(meshFile) 

global V mu0 nMatType nMatProp numNodeElem nMagnet lr_pbc bt_pbc DIM matProp

fid = fopen(meshFile);

meshType_str = fscanf(fid,'%s',1) ; 

% element types 

if strcmp(meshType_str, "HEX")
	elemType = 'HEX1' ; 
elseif strcmp(meshType_str , "TET") 
	elemType = 'TET1' ; 
elseif strcmp(meshType_str , "TET2") 
	elemType = 'TET2' ; 
elseif strcmp(meshType_str , "HEX2") 
	elemType = 'HEX2' ; 
else
	meshType_str 
	error("invalid mesh type");
end

% read material property
[fid,nMatType,nMatProp, matProp] =  readMatProp(fid);  

nMagnet =  nMatType - 1; 
% read nodal coordinates
[fid,nNode,xNode] = readNodalCoord(fid) ;


nElem = fscanf(fid,'%d',1);
elemOrder = fscanf(fid,'%d',1);

if elemOrder == 1
	if strcmp(elemType, 'HEX1')
		
		numNodeElem = 8 ; numNodeFace = 4 ; 
	
	elseif strcmp(elemType, 'TET1' )

		numNodeElem = 4 ; numNodeFace =	3 ; 

	else
		elemType
		error('not yet implemented') ; 
	end
elseif elemOrder == 2

	if strcmp(elemType, 'HEX2' )

		numNodeElem = 27 ; numNodeFace =	9 ; 

    elseif strcmp(elemType, 'TET2' )

		numNodeElem = 10 ; numNodeFace	=	6 ; 
    else
		elemType
		error('not yet implemented') ; 
	end
else
	elemOrder
	error('not yet implemented');
end

for p = 1:nElem
	indx(p,1:numNodeElem)	= fscanf(fid,'%d',numNodeElem);

	elemMatIndx(p)			= fscanf(fid,'%d',1);
end

pbcvec = fscanf(fid,'%d',3);

nBCsets = pbcvec(1) + pbcvec(2) + pbcvec(3); 

for p = 1:2*nBCsets

	numBCnodes(p) = fscanf(fid,'%d',1) ;

	for q = 1:numBCnodes(p)
		bcNode(p,q) = fscanf(fid,'%d',1);
	end
	% numBCtri(p) = fscanf(fid,'%d',1) ;
    % 
	% for q = 1:numBCtri(p)
	% 	dumy = fscanf(fid,'%d',numNodeFace) ;
	% 	dumy = fscanf(fid,'%d',numNodeElem);
	% end
end



pbcNode_lr(1:numBCnodes(1),1) = bcNode(1,1:numBCnodes(1)) ;
pbcNode_lr(1:numBCnodes(2),2) = bcNode(2,1:numBCnodes(2)) ;

pbcNode_bt(1:numBCnodes(3),1) = bcNode(3,1:numBCnodes(3)) ;
pbcNode_bt(1:numBCnodes(4),2) = bcNode(4,1:numBCnodes(4)) ;
 

[indx,pbcNode_lr,pbcNode_bt] = changeNodeIndexStart(nElem,indx,pbcNode_lr,pbcNode_bt);

for p = 1:size(pbcNode_lr,1)
    masterNode =  pbcNode_lr(p,1);
    trailNode = pbcNode_lr(p,2);
    % pbcNode_lr(:,3) = lr_pbc(1);
    % pbcNode_lr(:,4) = lr_pbc(2);
    % pbcNode_lr(:,5) = lr_pbc(3);
end

for p = 1:size(pbcNode_bt,1)
    masterNode =  pbcNode_bt(p,1);
    trailNode = pbcNode_bt(p,2);
    % pbcNode_bt(:,3) = bt_pbc(1);
    % pbcNode_bt(:,4) = bt_pbc(2);
    % pbcNode_bt(:,5) = bt_pbc(3);
end

pbcNode =  pbcNode_bt;
leftNode = pbcNode_lr(:,1);
rightNode =  pbcNode_lr(:,2);

% nPresDisp = size(pbcNode,1); 

%%%%%% preprocessor for magnets
% 
% %%%% magnet dimension
% mu0 = 4*pi*10^-7;
% d =  6e-3;
% l =  10e-3;
% V = pi*(d/2)^2*l;
% Br= 1.17; %%%% remenece is donated by Br
% m =  1/mu0*Br*V;
% 
% magnetInfo = cell(nMagnet,3);
% 
% elemMagnetIndx = elemMatIndx - 1;
% 
% for magnet = 1:nMagnet
%     CM_elem = [];
%     id = find(elemMagnetIndx ==  magnet);
%     magnetElem =  indx(id,1:numNodeElem);
%     CM_elem = cal_elem_CM(elemType,xNode,magnetElem,CM_elem);
%     magnetInfo{magnet,1} =  magnetElem;
%     magnetInfo{magnet,2} = CM_elem;
%     magnetInfo{magnet,3} = m;
%     m = -m;
% end
% 
% [targetElemID, cmList, etaList] = getCM_based_on_an_element(magnetInfo, indx, xNode, elemType);
% 
% magnetConnection = [1 2; 3 2; 3 4];

fclose(fid) ;

end



function [fid,nMatType,nMatProp,elemMatProp] =  readMatProp(fid) 

% material property: 
nMatType = fscanf(fid,'%d',1) ;
nMatProp = fscanf(fid,'%d',1) ; 

for p = 1:nMatType
	matIndx(p)		  = fscanf(fid,'%d',1);
	elemMatProp(p,1:nMatProp) = fscanf(fid,'%f',nMatProp) ;
end

end

function [fid,nNode,xNode] = readNodalCoord(fid)

nNode = fscanf(fid,'%d',1);
for p = 1:nNode
	xNode(p,1:3) = fscanf(fid,'%f',3);
end

end

function [indx,pbcNode_lr,pbcNode_bt] = changeNodeIndexStart(nElem,indx,pbcNode_lr,pbcNode_bt)

nNodeElem = size(indx,2);

for p = 1:nElem
	for q = 1:nNodeElem
		indx(p,q) = indx(p,q) + 1; 
	end
end
% 
% sz1 = size(pbcNode,1) ; 
% sz2 = size(pbcNode,2) ;
% 
% for p = 1:sz1
% 	for q = 1:2 
% 
% 		pbcNode(p,q) = pbcNode(p,q) +1 ; 
% 	end
% end

pbcNode_lr = pbcNode_lr + 1;
pbcNode_bt = pbcNode_bt + 1;

end

function magnetNode =  listMagnetNode(magnetElem,magnetNode)
for node = magnetElem(:)
    if ~ismember(node,magnetNode)
        magnetNode = [magnetNode;node];
    end
end
end

function CM_elem = cal_elem_CM(elemType,xNode,magnetElem,CM_elem)

global DIM
numNodeElem =  size(magnetElem,2);
[shapeFunction, nGaussPts, gaussPts, gaussWeights] = ...
getShapeFunction_and_gaussInfo(elemType);
for elem = 1:size(magnetElem,1)

	node  	=  	magnetElem(elem,:);
	Xel	=	xNode(node,:)'; 

%	CM = sum(xNode(node,:))./4;

	[CM,v] = getElemCM(DIM,Xel, shapeFunction, ...
			nGaussPts , gaussPts , gaussWeights) ;

	CM_elem = [CM_elem;[CM' v]];
end

end
