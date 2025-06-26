function [K_global,P_global]  = ...
    element_stiffness_assemble(nDim,elem_node, eq_num, ...
    elem_num,Kel,Pel, K_global,P_global)

% This function assembles the element stiffness matrices into global
% matricies, partitioned according to the equation numbering.


%%%%% force assembly
% count =  0;
% for nodeID =  1:size(elem_node,2)
%     node =  elem_node(elem_num,nodeID);
%     for dim = 1:nDim
%         count = count + 1;
%         P_global(nDim*(node-1)+dim) = P_global(nDim*(node-1)+dim) + Pel(count);
%     end
% end


for i = 1:size(elem_node,2) % nodes in each element
    node_i = elem_node(elem_num,i);
    for j = 1:nDim % 2 directions
        row = eq_num(node_i,j); % row in global system
        for k = 1:size(elem_node,2)  % nodes in each element
            node_k = elem_node(elem_num,k);
            for l = 1:nDim % 2 directions
                col = eq_num(node_k,l); % column in global system
                    K_global(row,col) = K_global(row,col) + Kel(nDim*(i-1)+j,nDim*(k-1)+l);
            end
        end
            P_global(row) = P_global(row) + Pel(nDim*(i-1)+j);
        
    end
end
 % end of function
