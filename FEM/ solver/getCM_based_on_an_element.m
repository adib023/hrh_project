function [targetElemID, cmList, etaList] = getCM_based_on_an_element(magnetInfo, indx, xNode, elemType)

global a1 a2

for magnet = 1: 4
    elemList = magnetInfo{magnet,1};
    XcmList = magnetInfo{magnet,2};
    vXcm = 0;
    v = 0;

    for elem = 1: size(elemList,1)
        xcm =  XcmList(elem,1:3);
        velem =  XcmList(elem,4);
        vXcm =  vXcm + velem * xcm;
        v = v + velem;
    end

    CM_actual(magnet,:) =  vXcm/v;
    V_actual(magnet,:) = v;
end

cm1 = (CM_actual(1,:)*V_actual(1) +  (CM_actual(3,:) - a1)*V_actual(3))/...
    (V_actual(1) + V_actual(3));
cm2 = (CM_actual(2,:)*V_actual(2) +  (CM_actual(4,:) - a2)*V_actual(4))/...
    (V_actual(2) + V_actual(4));

cmList = [cm1;cm2]

target_elemNode = [];
% iniCM = [];


for p = 1:2
    cm =  cmList(p,:);
    dist = 1;

    magElem = magnetInfo{p,1};
    magCM =   magnetInfo{p,2};

    for q = 1: size(magCM,1)

        target_cm =  magCM(q,1:3);

        if norm(cm -  target_cm) <= dist
            dist =  norm(cm - target_cm);
            id = q;
        end

    end
    target_elemNode = [target_elemNode; magElem(id,:)];
end


targetElemID  = [];

for p = 1 : 2
    elemNode = target_elemNode(p,:);
    for q = 1: size(indx,1)
        if sum(ismember(elemNode, indx(q,:))) == size(indx,2)
            targetElemID = [targetElemID;q];
            break
        end
    end
end



etaList = [];

for p = 1:2
    Xel = xNode(target_elemNode(p,:),:)';
    xcm =  cmList(p,:)';

    eta = getCM_refCoord(elemType,Xel,xcm);
    etaList = [etaList;eta'];
end




end




function [eta] = getCM_refCoord(elemType,Xel,xcm)

% inputs: 
%	shapeFunction: routine that evaluates shape function and derivatives
%	Xel: element nodal coordinates 
%	xcm: center of mass (point for which reference coord is needed) 

% outputs:
%	eta: reference coordinate

% set ref. coord to center of cell as initial guess
[shapeFunction, ~, ~, ~] = getShapeFunction_and_gaussInfo(elemType);

eta = zeros(3,1) ; 

err = 1.0 ; 
TOL = 1e-3 ; 
count = 0 ; max_iter = 20 ; 

while (err > TOL)

	[N,DN] = shapeFunction(eta) ; 

	% solve equation: X_i N_i(eta) = x_cm by Newton Raphson
    
	d_eta = (Xel*DN)\(xcm - Xel*N') ; 

	eta  = eta + d_eta ;

	if norm(d_eta) < TOL
		break ; 
	end

	count = count + 1 ; 
	if (count > max_iter)
		'no convergence'
		exit ;
	end
end

end




