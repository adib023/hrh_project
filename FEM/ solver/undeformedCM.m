%function [cm1,cm2,cm3,cm4] = undeformedCM(elemType,magnetInfo, xNode)

clc
clear all

global DIM V mu0 nMatType nMatProp numNodeElem nMagnet...
    a1 a2 lr_pbc bt_pbc

meshFile = '../meshFiles/UCmagnet_hex2_fem.txt';

DIM = 3;

a2 =  [0.02620 0 0];

a1 =  [cos(pi/3)*a2(1) sin(pi/3)*a2(1) 0] ; 

isSparse = 1 ; 
isMagnet = 1;

mu0 = 4*pi*10^-7;

lr_strain = [-0.1 0 0];
bt_strain = [0 -0.1 0];

lr_pbc =  a1.*lr_strain;
bt_pbc =  a2.*bt_strain;



[elemType,nNode, nElem, xNode, indx ,matProp,elemMatIndx,...
    pbcNode,magnetInfo,magnetConnection, elemMagnetIndx] = preProcess_new(meshFile); 
fixedNode = 18;

[nDof,nPBC,eqNum, matStiff,Wmat] =...
    initialize(nNode,pbcNode,matProp,fixedNode);



[elemType,nNode, nElem, xNode, indx ,matProp,elemMatIndx,...
    pbcNode,magnetInfo,magnetConnection, elemMagnetIndx] = preProcess_new(meshFile); 

[shapeFunction, nGaussPts, gaussPts, gaussWeights] = ...
getShapeFunction_and_gaussInfo(elemType);


for magnet = 1: 4
    elemList = magnetInfo{magnet,1};
    XcmList = magnetInfo{magnet,2};
    vXcm = 0;
    v = 0;

    for elem = 1: size(elemList,1)
        velem = 0;
        Xel = xNode(elemList(elem,:),:)';
        for p = 1:nGaussPts
            r =  gaussPts(p,1:3);
            w = gaussWeights(p);

            [N,DN] = shapeFunction(r);
            J = Xel*DN;
            detJ =  det(J);

            velem = velem +  detJ * w; 
        end
        v = v+ velem;
        vXcm =  vXcm + velem * XcmList(elem,:);
    end

    CM_actual(magnet,:) =  vXcm/v;
    V_actual(magnet,:) = v;
end

CM_actual
% V
% V_actual(1) + V_actual(3)
% V_actual(2) + V_actual(4)


% V_actual(3)
% V_actual(4)
cm1 = (CM_actual(1,:)*V_actual(1) +  (CM_actual(3,:) - a1)*V_actual(3))/...
    (V_actual(1) + V_actual(3))
cm2 = (CM_actual(2,:)*V_actual(2) +  (CM_actual(4,:) - a2)*V_actual(4))/...
    (V_actual(2) + V_actual(4))
% cm3 = (CM_actual(3,:)*V_actual(3) +  (CM_actual(1,:) + a1)*V_actual(1))/...
%     (V_actual(1) + V_actual(3))
% cm4 = (CM_actual(4,:)*V_actual(4) +  (CM_actual(2,:) + a2)*V_actual(2))/...
%     (V_actual(2) + V_actual(4))

% end



%%%%% find out which element belong to the center of mass

cm_list = magnetInfo{1,2};
elem_list =  magnetInfo{1,1};
dist = 1;
cm1 = cm1;


for p = 1:size(cm_list,1)
    cm = cm_list(p,:);
    if norm(cm-cm1) <= dist
        dist =  norm(cm-cm1);
        id = p;
        magElemNode =  elem_list(p,:);
    end
end
% dist
cm = cm_list(id,:)
cm1


%%%% plot element


figure

for p = 1:length(magElemNode)
    xc1 =  xNode(magElemNode(p),:);
    for q = p+1:length(magElemNode)
        xc2 = xNode(magElemNode(q),:);
        plot3([xc1(1) xc2(1)],[xc1(2) xc2(2)],[xc1(3) xc2(3)],'k-')
        hold on
    end
end
scatter3(cm1(1),cm1(2),cm1(3),100, "b*")
hold on
scatter3(cm(1),cm(2),cm(3),100, "r*")
hold on





% magElemNode
%%%%%% find the elemenets in indx

% targetElemList = [];
% for p = 1:nElem
%     if sum(ismember(indx(p,:),magElemNode)) > 0
%         targetElemList  =  [targetElemList p];
%     end
% end
% 
% targetElemList
% 
% 
% %%%%% cal cm for targeted element
% 
% v = 0;
% vxcm = 0;  
% for p = targetElemList
% 
%     node  	=  	indx(p,:);
% 	Xel	=	xNode(node,:)'; 
% 
%     [cm,vol] = getElemCM(DIM,Xel, shapeFunction, ...
% 			nGaussPts , gaussPts , gaussWeights) ;
% 
%     vxcm = vxcm +  cm * vol
%     v=  v + vol;
% end
% 
% xcm =  vxcm/v
% cm1

