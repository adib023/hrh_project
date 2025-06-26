function [UUR,normDisp,normForce] = solve(nDim,n_node, Wmat , K_global,P_global, UR_old)

%add_row = [];

Ured = K_global  \ P_global;  %%%%delta(u) = UF

% normDisp = norm(Ured) ; 
% normForce = norm(P_global) ; 

UF = Wmat * Ured ; 


% insert the free and prescribed displacements  into Uur
row =  0;
for i = 1:n_node
    for j = 1:nDim
        row = row + 1;
        UUR(i,j) = UR_old(i,j) - UF(row); 
    end
end

normDisp = norm(UF);
normForce =  norm(P_global);

return
