function [fid,strainList] =  NRsolver(isSparse, elemType,n_node, nElem ,...
    xNode, indx, elemMatIndx, matStiff,...
    eqNum , Wmat ,  nDim,...
    pbcNode,leftNode, rightNode,fid)

global TOL_DISP TOL_FORCE minDelta maxDelta maxItr

strainFid =  fopen('../Results/strainList.txt','w');

UUR = zeros(n_node,nDim);
UUR1 = UUR;
fid =  write(fid,UUR1);
TOL_DISP = 1e-6;
TOL_FORCE = 1e-6;

maxDelta =  0.01;
minDelta  = 1e-5;
delta = maxDelta; % starting delta
maxItr = 20;

dfctrCurrent = 0 ;
dfctr = maxDelta;
step = 0;
strainList = [];



while (dfctr <=  1)

    strcat('Step:  ',num2str(step+1),'  ', ...
        'detla: ',num2str(delta),'  ','dfctr: ',num2str(dfctr))

    UUR2 =  applyPBC(pbcNode,leftNode, rightNode, UUR1,delta);

    UUR2_ini = UUR2;
    converge= false;
    Energy1 = 0;
    itr =  0;

    while (converge == false)
        itr =  itr + 1

        [Energy, P_global, K_global] = ...
            total_assembly(isSparse , elemType,nElem, xNode, indx,...
            elemMatIndx, matStiff, eqNum, Wmat , UUR2, ...
            true,true);

        UR_old  = UUR2;
        if itr <= maxItr || itr>2*maxItr
            [UUR2,normDisp,normForce] = solve(nDim,n_node, Wmat , ...
        		K_global,P_global, UR_old);

            'change of energy:(-  decrease, + increase)'
            delEnergy =  Energy - Energy1
            Energy1 = Energy;
            [normDisp normForce]

        else

            %%%%% back tracking armijo search %%%%%%%
            [UUR2,normForce,normDisp, ispossible] = ...
                steepestDescent(isSparse,elemType,Wmat, ...
                UUR2_ini, nElem,...
                xNode, indx,elemMatIndx, matStiff, eqNum,normDisp);

            UUR2_ini = UUR2;
            %[normDisp normForce]
            if ~ispossible
                'Back Tracking not possible'
                itr = 3*maxItr + 1;
            end
        end
    	dfctrCurrent = dfctr ;
        [UUR1, converge,dfctr,delta,step,fid,writeStrain]  = ...
            makeDecision(normDisp,normForce,UUR1, n_node, ...
            K_global , Wmat , dfctr,delta,itr,step,fid,UUR2) ;

    end

    %write strain and energy to file at end of each converged step:
    strainList =  [step dfctrCurrent Energy]
    if writeStrain
        strainFid = write(strainFid,strainList);
    end
end

M_global = mass_matrix(elemType,nElem, ...
    xNode, indx, elemMatIndx, eqNum, UUR1, Wmat);


[~, ~, K_global] = ...
    total_assembly(isSparse , elemType,nElem, xNode, indx,...
    elemMatIndx, matStiff, eqNum, Wmat , UUR1, ...
    false,true);

[freq, modeshape] = determineMode(M_global, K_global, 5000, 200, Wmat, ...
    'defFreq.txt','defMode.txt');

alpha = 0.01;
writeMode(alpha , UUR1,freq,modeshape,'vtkDefMode.txt') 


fid_final_UUR =  fopen('../Results/final_UUR.txt','w');
fid_step = fopen('../Results/stepCount.txt','w');

fid_final_UUR = write(fid_final_UUR, UUR1);
fid_step = write(fid_step, strainList(1));

fclose(fid_step);
fclose(fid_final_UUR);
fclose(strainFid);

end


function UUR2  =  applyPBC(pbcNode,leftNode, rightNode, UUR1,delta)

%%% UUR1 is converged UUR from previous steo
%%% UUR2 is the UUR introduced with PBC to start new step of NRsolver

global left_disp right_disp

UUR2 = UUR1;
for p = 1:size(pbcNode,1)

    masterID =  pbcNode(p,1);
    followerID =  pbcNode(p,2);

    UUR2(followerID,1) = UUR2(masterID,1);
    UUR2(followerID,2) = UUR2(masterID,2);
    UUR2(followerID,3) = UUR2(masterID,3);
end

for p = 1: length(leftNode)
    UUR2(leftNode(p),:) = UUR2(leftNode(p),:) + delta * left_disp;
    UUR2(rightNode(p),:) = UUR2(rightNode(p),:) + delta * right_disp;
end

end

