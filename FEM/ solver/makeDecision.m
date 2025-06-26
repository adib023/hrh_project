function [UUR1, converge,dfctr,delta, step,fid, writeStrain]  = makeDecision(normDisp,normForce,UUR1, nNode ,  ...
				Kred , Wmat, dfctr,delta,itr,step,fid,UUR2)


%    nElem, xNode, indx,matStiff,...
%    pbcNode, fixID, eqNun, nDof, UUR2,converge,dfctr,delta,itr,step,fid, ...
%    magnetInfo, magnetConnection,elemMatIndx, elemMagnetIndx)

global TOL_DISP TOL_FORCE minDelta maxDelta maxItr

converge = false ; 
writeStrain = false;

if((normDisp < TOL_DISP) && (normForce < TOL_FORCE))
    converge = true;
            [minEval,minUF] = checkMinEval( Kred , Wmat );
            if minEval < 0
                writeStrain = false;
                "converged but unstable solution"
                [dfctr,delta] = controlDispFctr(dfctr,delta);
                if delta< minDelta
                    [UUR1,dfctr,delta] = perturb(nNode,UUR1,minUF,dfctr,delta);
                end
            else
                "converged"
                writeStrain = true;
                step = step + 1;
                fid = write(fid,UUR2);
                UUR1 = UUR2;
                remainDfctr = 1- dfctr;
                if remainDfctr < maxDelta && remainDfctr > minDelta
                    delta = remainDfctr;
                    dfctr = dfctr + delta;
                elseif remainDfctr < minDelta
                    dfctr  = 1.2;
                else
                    delta =  min([maxDelta 1.2*delta]);
                    dfctr = delta + dfctr;
                end
            end
end

if itr > 3*maxItr
    'Does not converge'
    [dfctr,delta] = controlDispFctr(dfctr,delta);
    if delta< minDelta
            [~,minUF] = checkMinEval( Kred , Wmat );
        [UUR1,dfctr,delta] = perturb(nNode,UUR1,minUF,dfctr,delta);%%%%% check this minUF
    end
    converge =  true; %%%%% it does not indicate the solution converged but used to teminate the while loop
    writeStrain = false;
end


end

function [UUR1,dfctr,delta] = perturb(nNode,UUR1,minUF,dfctr,delta)
perUUR =  zeros(nNode,3);
c = 0;
perdel = 10^-2 * norm(UUR1(:))/norm(minUF);
dfctr = dfctr - delta;
delta =  1e-4;
dfctr = dfctr + delta;
for node =  1:nNode
    for dim = 1:3
        c = c+1;
        perUUR(node,dim) = perdel*minUF(c);
    end
end

UUR1 = UUR1 + perUUR;
'Perturbed'

end

function writeModes()

modefile = fopen('mode','w') ; 
modefile = write(modefile,perUUR) ; 

end

function [minEval,minUF] = checkMinEval( Kred , Wmat )

[V,eigVals] =  eigs(Kred,5,0);

'The lowest eigenvalues'
eigVals =  diag(eigVals)
minEval = min(eigVals);

%nNode = size(xNode,1) ; 

% writing mode shapes for debugging
%writeModeShapes(V,eigVals , nDim , nNode , ...
%			eqNum,pbcNode,elim_eqn) ;


id =  find(minEval == eigVals);
minUF =  V(:,id);

minUF = Wmat * minUF ; 

end

function [fullMode] = getFullMode(mode,nDim,nNode, ... 
				eqNum,pbcNode,elim_eqn) 

T = eye(nDim*nNode);
T(:,elim_eqn)= [];

fullMode = T*mode ;

for p = 1:size(pbcNode,1)
    lead =  pbcNode(p,1);
    trail =  pbcNode(p,2);
    for q = 1:nDim
        l_eqn = eqNum(lead,q);
        t_eqn = eqNum(trail,q);
        fullMode(t_eqn) = fullMode(l_eqn);
    end
end


end


function writeModeShapes(V,eigVals , nDim , nNode , ...
			eqNum,pbcNode,elim_eqn)

ndof = nDim * nNode ; 

for p = 1:size(eigVals,1) 

	mode = V(:,p) ;  

	[fullMode] = getFullMode(mode,nDim,nNode, ... 
				eqNum,pbcNode,elim_eqn) ;  

	fid = fopen(strcat('mode_' , num2str(p)), 'w'); 

	for q = 1:ndof
		fprintf(fid,'%.16f \n',fullMode(q)) ;   
	end

	fclose(fid) ; 
end

messageOut = 'exiting in writeModeShapes'
exit ; 
end
