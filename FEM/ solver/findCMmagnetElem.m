function [targetElemID, iniCM] = findCMmagnetElem(magnetInfo,indx)

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

cmList = [cm1;cm2];
target_elemNode = [];
iniCM = [];
for p = 1:2
    cm =  cmList(p,:);
    dist = 1;

    magElem = magnetInfo{p,1};
    magCM =   magnetInfo{p,2};

    for q = 1: size(magCM,1)

        target_cm =  magCM(q,1:3);

        if norm(cm -  target_cm) <= dist
            iniCM =  [iniCM;target_cm];
            dist =  norm(cm - target_cm);
            elemNode =  magElem(q,:);
        end

    end
    target_elemNode = [target_elemNode; elemNode]; 
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

end


