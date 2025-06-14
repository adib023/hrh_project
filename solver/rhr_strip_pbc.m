function rhr_strip_pbc(N2,k,kt,pres_disp)
%%% re-ent lattice vector
a1r = [2*cos(pi/12) 0];
a2r = [cos(-pi/12) 1+sin(-pi/12)];
%%% hexa lattice vector
a1h = [2*cos(pi/12) 0];
a2h = 2*cos(pi/12)*[cos(pi/3) sin(pi/3)];

%%% initial coords for first re-ent
x1 = [0 0];
x2 = [cos(-pi/12) sin(-pi/12)];

numNode = 0;

for p = 0:N2-1
    node1 =  p*a2r + x1;
    node2 =  p*a2r + x2;
    numNode =  numNode + 2;
    ucNode(1,p+1,1:2) =  [numNode-1 numNode];
    xCoord(numNode-1,1:2) = node1;
    xCoord(numNode,1:2) = node2;
end

%%% initial coord for interface hexa
x1 = a2r + xCoord(ucNode(1,N2,1),1:2);
x2 =  x1 + cos(pi/12)/cos(pi/6)*[cos(pi/6) sin(pi/6)];

%%% define nodal coordinate for hexagonal lattice
node1 = a2h+x1;
node2 = a2h+x2;
numNode = numNode + 2;
ucNode(1,N2+1,1:2) = [numNode-1 numNode];
xCoord(numNode-1,1:2) = node1;
xCoord(numNode,1:2) = node2;

%%%% initial coord for re-ent2
x1 = xCoord(nodeIDperUC(1,N2+1,1),1:2) + a2h;
x2 = x1 +  [cos(-pi/12) sin(-pi/12)];

for p = 0:N2-1
    node1 =  p*a2r + x1;
    node2 =  p*a2r + x2;
    numNode =  numNode + 2;
    ucNode(1,N2+1+p+1,1:2) =  [numNode-1 numNode];
    xCoord(numNode-1,1:2) = node1;
    xCoord(numNode,1:2) = node2;
end

N2 =  N2 + N2 + 1;

numElem = 0;

for p = 1:N2
    numElem = numElem + 1;
    indx(numElem,1:3) = [nodeIDperUC(p,q,1) nodeIDperUC(p,q,2) k];
    if p ~= N2
        numElem = numElem + 1;
        indx(numElem,1:3) = [nodeIDperUC(p,q,2) nodeIDperUC(p,q+1,1) k];
    end
end

numRot = 0;
%%%% define rotational element

for p = 1:N2
    if p ~= 1
        numRot = numRot + 1;
        rotIndx(numRot,1:3) = [ucNode(1,p,1) ucNode(1,p-1,2) ucNode(1,p,2)];
    end
    if p ~= N2
        numRot = numRot + 1;
        rotIndx(numRot,1:3) = [ucNode(1,p,2) ucNode(1,p+1,1) ucNode(1,p,1)];
    end
end

rotIndx(:,4) = kt;

for p = 1: size(rotIndx,1)
    n1 = rotIndx(p,1);
    n2 = rotIndx(p,2);
    n3 = rotIndx(p,3);
    
    xc = xCoord([n1 n2 n3],:);
    
    costh_0 =  dot((xc(2,:)-xc(1,:)),(xc(3,:)-xc(1,:)))/...
        norm(xc(1,:)-xc(2,:))/norm(xc(1,:)-xc(3,:));
    sinth_0 =  norm(cross(([xc(2,:) 0]-[xc(1,:) 0]),([xc(3,:) 0]-[xc(1,:) 0])))/...
        norm(xc(1,:)-xc(2,:))/norm(xc(1,:)-xc(3,:));
    th_0 =  atan2(sinth_0,costh_0);
    rotIndx(p,5) = th_0;
end

%%%%% define periodic boundary condition
%%%%% in list first there will be right conditions
nAxPBC = 0; 
for p = 1:N2
    nAxPBC = nAxPBC + 1;
    n1 =  ucNode(1,p,2);
    n2 =  ucNode(1,p,1);
    l_0 = norm(xCoord(n1,:)-(xCoord(n2,:)+a1r));
    axPBC(nAxPBC,:) = [n1 n2 k [0 0] a1r l_0];
end
 
for p = 1:N2
    nAxPBC = nAxPBC + 1;
    n1 =  ucNode(1,p,1);
    n2 =  ucNode(1,p,2);
    l_0 = norm(xCoord(n1,:)-(xCoord(n2,:)-a1r));
    axPBC(nAxPBC,:) = [n1 n2 k [0 0] -a1r l_0];
end

%%%% define rotational pbc

nRotPBC = 0;

for p = 1:N2
    if p ~=N2
        nRotPBC = nRotPBC + 1;
        n1 = ucNode(1,p,2);
        n2 = ucNode(1,p+1,1);
        n3 = ucNode(1,p,1);
        l1 = [0 0];
        l2 = [0 0];
        l3 = a1r;
        rotPBC(numRotBC,:) = [n1 n2 n3 kt l1 l2 l3];
    end
    nRotPBC = nRotPBC + 1;
    n1 = ucNode(1,p,2);
    n2 = ucNode(1,p,1);
    n3 = ucNode(1,p,1);
    l1 = [0 0];
    l2 = [0 0];
    l3 = a1r;
    rotPBC(numRotBC,:) = [n1 n2 n3 kt l1 l2 l3];
    
    nRotPBC = nRotPBC + 1;
    n1 = ucNode(1,p,1);
    n2 = ucNode(1,p,2);
    n3 = ucNode(1,p,2);
    l1 = a1r;
    l2 = a1r;
    l3 = [0 0];
    rotPBC(numRotBC,:) = [n1 n2 n3 kt l1 l2 l3];
    
    if p ~= 1
        nRotPBC = nRotPBC + 1;
        n1 = ucNode(1,p,1);
        n2 = ucNode(1,p-1,2);
        n3 = ucNode(1,p,2);
        l1 = a1r;
        l2 = a1r;
        l3 = [0 0];
        rotPBC(numRotBC,:) = [n1 n2 n3 kt l1 l2 l3];
    end
    
end

for p = 1:N2
    if p ~=N2
        nRotPBC = nRotPBC + 1;
        n1 = ucNode(1,p,2);
        n2 = ucNode(1,p+1,1);
        n3 = ucNode(1,p,1);
        l1 = -a1r;
        l2 = -a1r;
        l3 = [0 0];
        rotPBC(numRotBC,:) = [n1 n2 n3 kt l1 l2 l3];
    end
    nRotPBC = nRotPBC + 1;
    n1 = ucNode(1,p,2);
    n2 = ucNode(1,p,1);
    n3 = ucNode(1,p,1);
    l1 = [0 0];
    l2 = [0 0];
    l3 = -a1r;
    rotPBC(numRotBC,:) = [n1 n2 n3 kt l1 l2 l3];
    
    nRotPBC = nRotPBC + 1;
    n1 = ucNode(1,p,1);
    n2 = ucNode(1,p,2);
    n3 = ucNode(1,p,2);
    l1 = -a1r;
    l2 = -a1r;
    l3 = [0 0];
    rotPBC(numRotBC,:) = [n1 n2 n3 kt l1 l2 l3];
    
    if p ~= 1
        nRotPBC = nRotPBC + 1;
        n1 = ucNode(1,p,1);
        n2 = ucNode(1,p-1,2);
        n3 = ucNode(1,p,2);
        l1 = a1r;
        l2 = a1r;
        l3 = [0 0];
        rotPBC(numRotBC,:) = [n1 n2 n3 kt l1 l2 l3];
    end
    
end


end