function freq_res_anal(kiri)

[K,f_int] = getStiffness_matrix(kiri.xCoord,kiri.indx,kiri.rotIndx,1,0);
[K,~] = apply_PBC(kiri.xCoord,kiri.axialPBC,kiri.rotPBC,K,f_int,1,0);

F =  zeros(2*kiri.nNode,1);

numPresID =  length(kiri.presDispID);

kiri.presDispID(numPresID/2+1:numPresID)
F(2*kiri.presDispID(numPresID/2+1:numPresID)) = 1;



ilimID  = [];

for p =  1: numPresID
    if kiri.presDispX(p) == 0 & kiri.presDispY(p) == 0
        ilimID = [ilimID 2*kiri.presDispID(p)-1 2*kiri.presDispID(p)];
    end
end

K(ilimID,:) = [];
K(:,ilimID) = [];
F(ilimID) = [];

M =  eye(size(K));

N1  = kiri.N1;
N2 =  kiri.N2;

% targetNode1 =  kiri.id1(1:N1,3);
% targetNode2 =  kiri.id1(1:N1,7);
% targetNode3 =  kiri.id1(1:N1,10);
% 
% layer1 =  [2*targetNode1-1 2*targetNode1];
% layer2 =  [2*targetNode2-1 2*targetNode2];
% layer3 =  [2*targetNode3-1 2*targetNode3];

targetNode1 =  [kiri.id1(1:N1,1:N2) kiri.id2(1:N1,1:N2)];
targetNode2 =  [kiri.id1(1:N1,7) kiri.id2(1:N1,7)];
targetNode3 =  [kiri.id1(1:N1,10) kiri.id2(1:N1,10)];

layer1 =  [2*targetNode1-1 2*targetNode1];
layer2 =  [2*targetNode2-1 2*targetNode2];
layer3 =  [2*targetNode3-1 2*targetNode3];



datasheet = [];

for w =  0:0.0001:3.5
    w
    U =  (-w^2*M+K)^-1*F;
    T = eye(2*kiri.nNode);
    T(:,ilimID) = [];
    U = T*U;
    
    normU1 = sqrt(sum(U(layer1(:)).*conj(U(layer1(:)))))/norm(F);
    normU2 = sqrt(sum(U(layer2(:)).*conj(U(layer2(:)))))/norm(F);
    normU3 = sqrt(sum(U(layer3(:)).*conj(U(layer3(:)))))/norm(F);
    
    result = [w normU1 normU2 normU3];
    
    datasheet = [datasheet; result]; 
    
end 


file = fopen('../result/freq_res_hexa.txt','w');

write('../result/freq_res_hexa.txt',datasheet);

fclose(file);

end