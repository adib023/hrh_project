function [UUR2,normForce,normDisp, ispossible] = steepestDescent(isSparse,elemType, Wmat, UUR2_ini,nElem,...
    xNode, indx,elemMatIndx,matStiff, eqNum, normDisp)


UUR2= UUR2_ini;

[Energy1, Force, ~] = ...
            total_assembly(isSparse, elemType,nElem, xNode, indx,...
            elemMatIndx, matStiff, eqNum, Wmat, UUR2, true,false);


'energy for initial guess at the starting of steepest descent'
Energy1
norm0 = norm(Force)

Force = Wmat * Force;
ratio0 = 1;

for x =0:20
    x
    UUR3 =  UUR2;
    alpha =  0.1^x;
    deltaF = -alpha*Force;
    c = 0;
    for node = 1:size(xNode,1)
        for dim =  1: 3
            c = c+1;
            UUR3(node,dim) = UUR3(node,dim) + deltaF(c);
        end
    end

    [Energy, Force1, ~] = ...
            total_assembly(isSparse, elemType,nElem, xNode, indx,...
            elemMatIndx, matStiff,eqNum, Wmat, UUR3, true, false);


    norm1 = norm(Force1);
    ratio = norm1/norm0
        'Change in energy:'
        delEnergy =  Energy - Energy1
        Energy1= Energy;
    if norm1 < 0.9*norm0 
        UUR2 = UUR3;
        normForce = norm1;
        normDisp = norm(deltaF);
        ispossible =  true;
        break
    elseif ratio <= 1
	if ratio > ratio0
		UUR2 = UUR3_old;
		normForce = ratio0*norm0;
		normDisp = norm(deltaF_old);
		ispossible = true;
		break
	else
		ratio0 = ratio;
		UUR3_old = UUR3;
		deltaF_old = deltaF;
		ispossible = false;
		normForce =  norm0;
	end
    else
%         UUR2 =  UUR2;
	UUR3_old =  UUR3;
	deltaF_old = deltaF;
        ispossible = false;
        normForce =  norm0;
        %normDisp = normDisp;
    end
end
end
