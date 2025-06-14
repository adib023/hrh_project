function K_disp = dispK(xCoord,nodeIDperUC,K,vecK)

%vecK = [kx ky];

%[xCoord, nodeIDperUC, edgeIndx, rotIndx] =  dispAnalLattice(ka,kt);

%K =  getStiffness_matrix(xCoord,edgeIndx,rotIndx);

%%%%%% nodes at (0,0) unit cell%%%%%%%%

n1_0_0 = nodeIDperUC(2,2,1);
n2_0_0 = nodeIDperUC(2,2,2);

K = K(2*n1_0_0-1:2*n2_0_0,:);



K_disp = zeros(4);

for p = 1:3
    for q = 1:3
        n1 =  nodeIDperUC(p,q,1);
        n2 = nodeIDperUC(p,q,2);
        a1 = xCoord(n1,:) - xCoord(n1_0_0,:);
        a2 = xCoord(n2,:) - xCoord(n2_0_0,:);
        
        K(:,2*n1-1:2*n1) = K(:,2*n1-1:2*n1)*exp(i*dot(vecK,a1));
        K(:,2*n2-1:2*n2) = K(:,2*n2-1:2*n2)*exp(i*dot(vecK,a2));
        
        K_disp = K_disp + K(:,2*n1-1:2*n2);
    end
end


end

