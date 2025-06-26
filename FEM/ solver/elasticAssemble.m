function [TotalEnergy,P_global,K_global] = ...
    elasticAssemble (shapeFunction, ...
    nGaussPts, gaussPts, gaussWeights,...
    n_elem, coords, indx, mat_stiff,elemMatIndx, eq_num, UUR, ...
    TotalEnergy, P_global, K_global, calP, calK)

global DIM;

ndim = DIM;

% isSparse =1


fcn = @plus ;

nDof = size(K_global,1);


% P_g = zeros(nDof,1) ;
n_elem
parfor elem_num =  1:n_elem

    matID =  elemMatIndx(elem_num);

    % if matID  == 1
        elem_stiff = mat_stiff(1,:);

    % 
    % else
    %     elem_stiff = mat_stiff(2,:);
        
    % end
    elem_node = indx(elem_num,:);

    [Kel, Pel,ElemEnergy,detVol]   =  element(shapeFunction, ...
        nGaussPts, gaussPts, gaussWeights,...
        ndim, elem_num, coords, indx, ...
        elem_stiff, [], UUR, calP, calK, 0);

    % if (detVol <= 0)
    % 	exit;
    % end

    [rowIndxVec,colIndxVec,fIndxVec] = elementAssembly(0,elem_node, ...
        eq_num) ;

    if calK
        K_global = fcn(K_global , sparse(rowIndxVec,colIndxVec,Kel,nDof,nDof));
    end
    if calP
        P_global	= fcn(P_global , sparse(fIndxVec,ones(size(fIndxVec)),Pel,nDof,1));
    end
    %	Pg(rowIndxVec) = Pg(rowIndxVec) + Pel ;

    TotalEnergy = TotalEnergy + ElemEnergy ;

    % if calK
    % K_global = Wmat'*K_global*Wmat ;   %%%%%% reduced K_global
    % end
    %
    %
    % if calP
    % P_global = Wmat'*P_global ; %%%%% reduced P_global
    % end

end
end