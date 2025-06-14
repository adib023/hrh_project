function [K,f] = getStiffness_matrix (xCoord,indx,rotIndx,flag1,flag2)

%%%% flag1 for stiffness matrix
%%%% flag2 for force vector

nNode =  size(xCoord,1);
nDofNode = 2;
nDofTotal =  nNode * nDofNode;

K = zeros(nDofTotal);
f = zeros(nDofTotal,1);

%%%%%% axial spring assembly%%%%%%%%%%
nElem = size(indx,1);
for elem = 1:nElem
    n1 = indx(elem,1);
    n2 = indx(elem,2);
    xc = xCoord([n1,n2],:);
    k = indx(elem,3);
    l_0 = indx(elem,4);
    
    
    if flag1 ==  1
        Kspr = getAxialSpringK(xc,k,l_0);
    else
        Kspr = 0;
    end
    
    if flag2 == 1
        fspr =  getAxialForce(xc,k,l_0);
    else
        fspr = 0;
    end
        
    [K,f] =  assembleLinearSprings(K,Kspr,f,fspr,n1,n2,flag1,flag2);
end

%%%%%%% roattional spring%%%%%%%%%%%
numR =  size(rotIndx,1);

for n =  1:numR
    n1 = rotIndx(n,1);
    n2 = rotIndx(n,2);
    n3 = rotIndx(n,3);
    
    xc(1,1:2) = xCoord(n1,1:2);
    xc(2,1:2) = xCoord(n2,1:2);
    xc(3,1:2) = xCoord(n3,1:2);
    
    k =  rotIndx(n,4);
    th_0 = rotIndx(n,5);
    eta = cross(([xc(2,:) 0]-[xc(1,:) 0]),([xc(3,:) 0]-[xc(1,:) 0]))/...
        norm(xc(1,:)-xc(2,:))/norm(xc(1,:)-xc(3,:));
    
    if eta(3) < 0
        temp = xc(2,:);
        xc(2,:) = xc(3,:);
        xc(3,:) = temp;
        temp_n =  n2;
        n2 = n3;
        n3 = temp_n;
    end
    
    if flag1 == 1
        Kspr = getRotK(xc,k,th_0);
    else
        Kspr = 0;
    end
    
    if flag2 == 1
        fspr = getTorForce(xc,k,th_0);
    else
        fspr = 0;
    end
    
    [K,f] =  assembleRotationalSprings(K,Kspr,f,fspr,n1,n2,n3,flag1,flag2);
    
end

end

function [K,f] =  assembleLinearSprings(K,Kspr,f,fspr,node1,node2,flag1,flag2)
if flag1 == 1
K(2*node1-1:2*node1,2*node1-1:2*node1) = K(2*node1-1:2*node1,2*node1-1:2*node1) + Kspr(1:2,1:2); 
K(2*node1-1:2*node1,2*node2-1:2*node2) = K(2*node1-1:2*node1,2*node2-1:2*node2) + Kspr(1:2,3:4);
K(2*node2-1:2*node2,2*node1-1:2*node1) = K(2*node2-1:2*node2,2*node1-1:2*node1) + Kspr(3:4,1:2);
K(2*node2-1:2*node2,2*node2-1:2*node2) = K(2*node2-1:2*node2,2*node2-1:2*node2) + Kspr(3:4,3:4);
end

if flag2 ==1
f(2*node1-1,:) = f(2*node1-1,:) + fspr(1);
f(2*node1,:) = f(2*node1,:) + fspr(2);
f(2*node2-1,:) = f(2*node2-1,:) + fspr(3);
f(2*node2,:) = f(2*node2,:) + fspr(4);
end
end

function [K,f] =  assembleRotationalSprings(K,kspr,f,fspr,n1,n2,n3,flag1,flag2)
if flag1 == 1
K(2*n1-1:2*n1,2*n1-1:2*n1) = K(2*n1-1:2*n1,2*n1-1:2*n1) + kspr(1:2,1:2);
K(2*n1-1:2*n1,2*n2-1:2*n2) = K(2*n1-1:2*n1,2*n2-1:2*n2) + kspr(1:2,3:4);
K(2*n1-1:2*n1,2*n3-1:2*n3) = K(2*n1-1:2*n1,2*n3-1:2*n3) + kspr(1:2,5:6);

K(2*n2-1:2*n2,2*n1-1:2*n1) = K(2*n2-1:2*n2,2*n1-1:2*n1) + kspr(3:4,1:2);
K(2*n2-1:2*n2,2*n2-1:2*n2) = K(2*n2-1:2*n2,2*n2-1:2*n2) + kspr(3:4,3:4);
K(2*n2-1:2*n2,2*n3-1:2*n3) = K(2*n2-1:2*n2,2*n3-1:2*n3) + kspr(3:4,5:6);

K(2*n3-1:2*n3,2*n1-1:2*n1) = K(2*n3-1:2*n3,2*n1-1:2*n1) + kspr(5:6,1:2);
K(2*n3-1:2*n3,2*n2-1:2*n2) = K(2*n3-1:2*n3,2*n2-1:2*n2) + kspr(5:6,3:4);
K(2*n3-1:2*n3,2*n3-1:2*n3) = K(2*n3-1:2*n3,2*n3-1:2*n3) + kspr(5:6,5:6);
end

if flag2 == 1
f(2*n1-1,:) = f(2*n1-1,:) + fspr(1);
f(2*n1,:) = f(2*n1,:) + fspr(2);
f(2*n2-1,:) = f(2*n2-1,:) + fspr(3);
f(2*n2,:) = f(2*n2,:) + fspr(4);
f(2*n3-1,:) = f(2*n3-1,:) + fspr(5);
f(2*n3,:) = f(2*n3,:) + fspr(6);
end
end
