function [Ng1, Ng2] = get_global_NoI_cm(shapeFunction,etaList, targetElemID, indx, xNode)

global DIM

[N_cm1, ~] =  shapeFunction(etaList(1,:));
NoI_cm1 =  kron(N_cm1,eye(DIM));

[N_cm2, ~] =  shapeFunction(etaList(2,:));
NoI_cm2 =  kron(N_cm2,eye(DIM));

Ng1 = zeros(DIM, DIM*size(xNode,1));
Ng2 = zeros(DIM, DIM*size(xNode,1));


elemNode = indx(targetElemID(1),:);
Ng1 =  assembly(elemNode, NoI_cm1 , Ng1);


elemNode = indx(targetElemID(2),:);
Ng2 =  assembly(elemNode, NoI_cm2 , Ng2);

Ng1 = sparse(Ng1);
Ng2 = sparse(Ng2);

end




function N_global = assembly(elemNode, N, N_global)

global DIM

numElemNode  =  length(elemNode);

for n =  1:numElemNode
   node = elemNode(n);
   N_global(:,DIM*node-DIM + 1:DIM*node) = N_global(:,DIM*node-DIM + 1:DIM*node)... 
                                    + N(:,3*n- DIM + 1:3*n);
end

end