function M_global = mass_matrix(elemType,nElem, ...
    xNode, indx, elemMatIndx, eqNum, UUR, Wmat)

global DIM matProp;
% matProp
ndof =  DIM * size(xNode, 1);

[shapeFunction, nGaussPts, gaussPts, gaussWeights] = ...
    getShapeFunction_and_gaussInfo(elemType);

fcn = @plus;

M_global =  sparse(ndof, ndof);

parfor elem =  1: nElem
    % elem
    elemNode =  indx(elem,:);
   % 
   % if elemMatIndx(elem) == 1
        rho =  matProp(1,3);
   % else
   %     rho = matProp(2,3);
   % end
    
  MEL =  element(shapeFunction, ...
        nGaussPts, gaussPts, gaussWeights,...
        DIM , elem, xNode, indx, ...
        [], rho, UUR, 0, 0, 1);
        

  [rowIndxVec,colIndxVec,~] = elementAssembly(0,elemNode, ...
        eqNum) ;

  M_global = fcn(M_global , sparse(rowIndxVec,colIndxVec,MEL,ndof,ndof));

end


M_global =  Wmat' * M_global * Wmat;

end
