function disp_anal(ka,kt,plotType,lattice)

a1 = [2*cos(pi/12) 0];
if lattice == 're'
a2 = [cos(-pi/12) 1+sin(-pi/12)]; %%% for re-entrat
elseif lattice == 'he'
a2 = 2*cos(pi/12)*[cos(pi/3) sin(pi/3)]; %%% for hexa
end

x1 = [0 0];
if lattice == 're'
x2 = [cos(-pi/12) sin(-pi/12)]; %%% for re-entrat
elseif lattice == 'he'
x2 = cos(pi/12)/cos(pi/6)*[cos(pi/6) sin(pi/6)]; %%% for hexa
end

[xCoord, nodeIDperUC, edgeIndx, rotIndx] =  dispAnalLattice(a1,a2,x1,x2,ka,kt);
K =  getStiffness_matrix(xCoord,edgeIndx,rotIndx,1,0);
% 
if plotType == "IBZ"
[OA,AB,BC,CD,DO]=getIBZcoord(a1,a2);

vecK_list = [OA AB BC CD DO];


for n = 1:size(vecK_list,2)
    K_disp = dispK(xCoord,nodeIDperUC,K,vecK_list(:,n));
    w2 = eig(K_disp);
    w = sqrt(w2);
    w = real(w);
    w = sort(w);
    wList(:,n) = w;
end


figure
for n = 1: size(vecK_list,2)
    wSize = length(wList(:,n));
    plot(n*ones(wSize,1),wList(:,n),"k.","markersize",5)
    hold on
end
xlim([1 n]);

print('-dpng','-r300',"../result/disp_anal_IBZ")


elseif plotType == "surf" || plotType == "cont"

%%%%%%%%%%%%%%% surface ploting from here%%%%%%%%%%%%%%

kx =  linspace(-2*pi,2*pi,201);
ky =  linspace(-2*pi,2*pi,201);

[Kx,Ky] = meshgrid(kx,ky);

for p = 1:size(Kx,2)
    for q = 1:size(Ky,1)
        vecK = [Kx(q,p),Ky(q,p)];
        K_disp = dispK(xCoord,nodeIDperUC,K,vecK);
        w2 = eig(K_disp);
        w = sqrt(w2);
        w = real(w);
        w = sort(w);
        W1(q,p) = w(1);
        W2(q,p) = w(2);
        W3(q,p) = w(3);
        W4(q,p) = w(4);
    end
end
% 
% 
%
if plotType == "surf"
    figure
    colormap jet
    surf(Kx,Ky,W1,"EdgeColor","none");
    hold on
%    colorbar
%     view(0,90)
%     title("First Band")
%     print('-dpng','-r300',"../result/firstband")
    %figure
    %colormap jet
    surf(Kx,Ky,W2,"EdgeColor","none");
    hold on
    %colorbar
    %view(0,90)
    %title("Second Band")
    %print('-dpng','-r300',"../result/secondband")
    %figure
    %colormap jet
    surf(Kx,Ky,W3,"EdgeColor","none");
    hold on
%     colorbar
%     view(0,90)
%     title("Third Band")
%     print('-dpng','-r300',"../result/thirdband")
%     figure
%     colormap jet
    surf(Kx,Ky,W4,"EdgeColor","none");
    hold on
%     colorbar
%     view(0,90)
%     title("Fourth Band")
%     print('-dpng','-r300',"../result/fourthband")
elseif plotType == "cont"
    fsize =  20;
    figure
    colormap jet
    contour(Kx,Ky,W1);
    hold on
    contour(Kx,Ky,W1,[0.26879, 0.26879], 'k--', 'LineWidth',2);
    set(gca,'fontsize',fsize)
    set(gca,'fontname','Times New Roman')
    %title("First Band")
    colorbar
    xlabel("$k_x$","interpreter","latex",'fontsize',fsize)
    ylabel("$k_y$","interpreter","latex",'fontsize',fsize)
    print('-depsc','-r300',"../result/band1")
    print('-dpng','-r300',"../result/band1")
    
    figure
    colormap jet
    contour(Kx,Ky,W2);
    set(gca,'fontsize',fsize)
    set(gca,'fontname','Times')
    title("Second Band")
    colorbar
    xlabel("$k_x$","interpreter","latex",'fontsize',fsize)
    ylabel("$k_y$","interpreter","latex",'fontsize',fsize)
    print('-depsc','-r300',"../result/band2")
    print('-dpng','-r300',"../result/band2")
    % 
    figure
    colormap jet
    contour(Kx,Ky,W3);
    set(gca,'fontsize',fsize)
    set(gca,'fontname','Times')
    title("Third Band")
    colorbar
    xlabel("$k_x$","interpreter","latex",'fontsize',fsize)
    ylabel("$k_y$","interpreter","latex",'fontsize',fsize)
    print('-depsc','-r300',"../result/band3")
    print('-dpng','-r300',"../result/band3")

    figure
    colormap jet
    contour(Kx,Ky,W4);
    set(gca,'fontsize',fsize)
    set(gca,'fontname','Times')
    title("Fourth Band")
    colorbar
    xlabel("$k_x$","interpreter","latex",'fontsize',fsize)
    ylabel("$k_y$","interpreter","latex",'fontsize',fsize)
    print('-depsc','-r300',"../result/band4")
    print('-dpng','-r300',"../result/band4")
end
end
 





end

function[coord, nodeIDperUC, edgeIndx, rotIndx] =  dispAnalLattice(a1,a2,x1,x2,k,k_rot)
%%%%% this routine generates 3by3 hexagonal lattice file for dispersion analysis%%%%
N1 = 3;
N2 = 3;
% 
% %%%%%lattice Vector%%%%%%
% a1 = [1 0];
% a2 = [cos(pi/3) sin(pi/3)];
% %%%%% lattice vectors a1, a2
% 
% %%%% local vector of two nodes in a unit cell
% x1 = 1/sqrt(3) * [cos(-pi/6) sin(-pi/6)];
% x2 = 2*x1;



%%%%%%%% This section creates a coordinate lists and 
%%%%%%%% assign nodeID (p+1,q+1)th unit cell in nodeIDperUC 3D matrix
nodeCount =  0;

for p  =  0:N1-1
    for q = 0:N2-1
        node1 = p * a1 + q*a2 + x1;
        node2 = p * a1 + q*a2 + x2;
        nodeCount =  nodeCount + 2;
        nodeIDperUC(p+1,q+1,1:2) = [nodeCount-1 nodeCount];
        coord(nodeCount-1,1:2) = node1;
        coord(nodeCount,1:2) = node2;
    end
end

numNode =  nodeCount;

%%%%%%% This section edge list that contains nodeID of each edge%%%%%%%
edgeCount = 0;
for p = 1:N1-1
    for q = 1:N2-1
        edgeCount =  edgeCount + 3;
        edgeIndx(edgeCount-2,1:2) = [nodeIDperUC(p,q,1) nodeIDperUC(p,q,2)]; 
        edgeIndx(edgeCount-1,1:2) = [nodeIDperUC(p,q,2) nodeIDperUC(p+1,q,1)];
        edgeIndx(edgeCount,1:2) =  [nodeIDperUC(p,q,2) nodeIDperUC(p,q+1,1)];
    end
end

for p = 1:N1-1
    for q = N2
        edgeCount =  edgeCount + 2;
        edgeIndx(edgeCount-1,1:2) = [nodeIDperUC(p,q,1) nodeIDperUC(p,q,2)]; 
        edgeIndx(edgeCount,1:2) = [nodeIDperUC(p,q,2) nodeIDperUC(p+1,q,1)];
    end
end

for p = N1
    for q = 1:N2-1
        edgeCount =  edgeCount + 2;
        edgeIndx(edgeCount-1,1:2) = [nodeIDperUC(p,q,1) nodeIDperUC(p,q,2)]; 
        edgeIndx(edgeCount,1:2) = [nodeIDperUC(p,q,2) nodeIDperUC(p,q+1,1)];
    end
end

edgeCount =  edgeCount + 1 ;
edgeIndx(edgeCount,1:2) = [nodeIDperUC(N1,N2,1) nodeIDperUC(N1,N2,2)];
edgeIndx(:,3) = k;

numElem =  edgeCount;

for p = 1:numElem
    n1 =  edgeIndx(p,1);
    n2 =  edgeIndx(p,2);
    edgeIndx(p,4) = norm(coord(n1,:)-coord(n2,:));
end

%%%% generate torsional %%%%%

numRot =  0;

for p = 2:N1-1
    for q = 2:N2-1
        numRot =  numRot + 1;
        rotIndx(numRot,1:3) = [nodeIDperUC(p,q,1)...
                               nodeIDperUC(p,q-1,2)...
                               nodeIDperUC(p,q,2)];
        numRot =  numRot + 1;
        rotIndx(numRot,1:3) = [nodeIDperUC(p,q,1)...
                               nodeIDperUC(p,q,2)...
                               nodeIDperUC(p-1,q,2)];
        numRot =  numRot + 1;                   
        rotIndx(numRot,1:3) = [nodeIDperUC(p,q,1)...
                               nodeIDperUC(p-1,q,2)...
                               nodeIDperUC(p,q-1,2)];
        numRot =  numRot + 1;
        rotIndx(numRot,1:3) = [nodeIDperUC(p,q,2)...
                               nodeIDperUC(p,q,1)...
                               nodeIDperUC(p+1,q,1)];
        numRot =  numRot + 1;                   
        rotIndx(numRot,1:3) = [nodeIDperUC(p,q,2)...
                               nodeIDperUC(p,q+1,1)...
                               nodeIDperUC(p,q,1)];
        numRot =  numRot + 1;                   
        rotIndx(numRot,1:3) = [nodeIDperUC(p,q,2)...
                               nodeIDperUC(p+1,q,1)...
                               nodeIDperUC(p,q+1,1)];
    end
end

for p = 1
    for q = 2:N2-1 
    numRot =  numRot + 1;
    rotIndx(numRot,1:3) = [nodeIDperUC(p,q,1)...
                               nodeIDperUC(p,q,2)...
                               nodeIDperUC(p,q-1,2)];
   numRot =  numRot + 1;
        rotIndx(numRot,1:3) = [nodeIDperUC(p,q,2)...
                               nodeIDperUC(p,q,1)...
                               nodeIDperUC(p+1,q,1)];
        numRot =  numRot + 1;                   
        rotIndx(numRot,1:3) = [nodeIDperUC(p,q,2)...
                               nodeIDperUC(p,q+1,1)...
                               nodeIDperUC(p,q,1)];
        numRot =  numRot + 1;                   
        rotIndx(numRot,1:3) = [nodeIDperUC(p,q,2)...
                               nodeIDperUC(p+1,q,1)...
                               nodeIDperUC(p,q+1,1)];
    end
end 

for p = 2:N1-1
    for q =1 
    numRot =  numRot + 1;
    rotIndx(numRot,1:3) = [nodeIDperUC(p,q,1)...
                               nodeIDperUC(p-1,q,2)...
                               nodeIDperUC(p,q,2)];
   numRot =  numRot + 1;
        rotIndx(numRot,1:3) = [nodeIDperUC(p,q,2)...
                               nodeIDperUC(p,q,1)...
                               nodeIDperUC(p+1,q,1)];
        numRot =  numRot + 1;                   
        rotIndx(numRot,1:3) = [nodeIDperUC(p,q,2)...
                               nodeIDperUC(p,q+1,1)...
                               nodeIDperUC(p,q,1)];
        numRot =  numRot + 1;                   
        rotIndx(numRot,1:3) = [nodeIDperUC(p,q,2)...
                               nodeIDperUC(p+1,q,1)...
                               nodeIDperUC(p,q+1,1)];
    end
end 

for p = N1
    for q = 2:N2-1 
     numRot =  numRot + 1;                   
        rotIndx(numRot,1:3) = [nodeIDperUC(p,q,2)...
                               nodeIDperUC(p,q+1,1)...
                               nodeIDperUC(p,q,1)];
    numRot =  numRot + 1;
        rotIndx(numRot,1:3) = [nodeIDperUC(p,q,1)...
                               nodeIDperUC(p,q-1,2)...
                               nodeIDperUC(p,q,2)];
        numRot =  numRot + 1;
        rotIndx(numRot,1:3) = [nodeIDperUC(p,q,1)...
                               nodeIDperUC(p,q,2)...
                               nodeIDperUC(p-1,q,2)];
        numRot =  numRot + 1;                   
        rotIndx(numRot,1:3) = [nodeIDperUC(p,q,1)...
                               nodeIDperUC(p-1,q,2)...
                               nodeIDperUC(p,q-1,2)];
    end
end 

for p = 2:N1-1
    for q =N2
   numRot =  numRot + 1;
        rotIndx(numRot,1:3) = [nodeIDperUC(p,q,2)...
                               nodeIDperUC(p,q,1)...
                               nodeIDperUC(p+1,q,1)];
    numRot =  numRot + 1;
        rotIndx(numRot,1:3) = [nodeIDperUC(p,q,1)...
                               nodeIDperUC(p,q-1,2)...
                               nodeIDperUC(p,q,2)];
        numRot =  numRot + 1;
        rotIndx(numRot,1:3) = [nodeIDperUC(p,q,1)...
                               nodeIDperUC(p,q,2)...
                               nodeIDperUC(p-1,q,2)];
        numRot =  numRot + 1;                   
        rotIndx(numRot,1:3) = [nodeIDperUC(p,q,1)...
                               nodeIDperUC(p-1,q,2)...
                               nodeIDperUC(p,q-1,2)];
    end
end 

numRot =  numRot + 1;
rotIndx(numRot,1:3) = [nodeIDperUC(1,1,2)...
    nodeIDperUC(1,1,1)...
    nodeIDperUC(2,1,1)];
numRot =  numRot + 1;
rotIndx(numRot,1:3) = [nodeIDperUC(1,1,2)...
    nodeIDperUC(1,2,1)...
    nodeIDperUC(1,1,1)];
numRot =  numRot + 1;
rotIndx(numRot,1:3) = [nodeIDperUC(1,1,2)...
    nodeIDperUC(2,1,1)...
    nodeIDperUC(1,2,1)];

numRot =  numRot + 1;
rotIndx(numRot,1:3) = [nodeIDperUC(N1,N2,1)...
    nodeIDperUC(N1,N2-1,2)...
    nodeIDperUC(N1,N2,2)];
numRot =  numRot + 1;
rotIndx(numRot,1:3) = [nodeIDperUC(N1,N2,1)...
    nodeIDperUC(N1,N2,2)...
    nodeIDperUC(N1-1,N2,2)];
numRot =  numRot + 1;
rotIndx(numRot,1:3) = [nodeIDperUC(N1,N2,1)...
    nodeIDperUC(N1-1,N2,2)...
    nodeIDperUC(N1,N2-1,2)];

numRot =  numRot + 1;
rotIndx(numRot,1:3) = [nodeIDperUC(1,N2,1)...
    nodeIDperUC(1,N2-1,2)...
    nodeIDperUC(1,N2,2)];
numRot =  numRot + 1;
rotIndx(numRot,1:3) = [nodeIDperUC(1,N2,2)...
    nodeIDperUC(1,N2,1)...
    nodeIDperUC(2,N2,1)];

numRot =  numRot + 1;
rotIndx(numRot,1:3) = [nodeIDperUC(N1,1,1)...
    nodeIDperUC(N1,1,2)...
    nodeIDperUC(N1-1,1,2)];
numRot =  numRot + 1;
rotIndx(numRot,1:3) = [nodeIDperUC(N1,1,2)...
    nodeIDperUC(N1,2,1)...
    nodeIDperUC(N1,1,1)];


rotIndx(:,4) = k_rot;

for elem = 1:size(rotIndx,1)
    n =  rotIndx(elem,1:3);
    xc =  coord(n,:);
    eta = cross(([xc(2,:) 0]-[xc(1,:) 0]),([xc(3,:) 0]-[xc(1,:) 0]))/...
        norm(xc(1,:)-xc(2,:))/norm(xc(1,:)-xc(3,:));
    
    if eta(3) < 0
%         eta
        temp = rotIndx(elem,2);
        rotIndx(elem,2) = rotIndx(elem,3);
        rotIndx(elem,3) = temp;
%         n =  rotIndx(elem,1:3);
%         xc =  xCoord(n,:);        
%         eta = cross(([xc(2,:) 0]-[xc(1,:) 0]),([xc(3,:) 0]-[xc(1,:) 0]))/...
%             norm(xc(1,:)-xc(2,:))/norm(xc(1,:)-xc(3,:))
        
    end
end


for p = 1: size(rotIndx,1)
    n1 = rotIndx(p,1);
    n2 = rotIndx(p,2);
    n3 = rotIndx(p,3);
    
    xc = coord([n1 n2 n3],:);
    
    costh_0 =  dot((xc(2,:)-xc(1,:)),(xc(3,:)-xc(1,:)))/...
        norm(xc(1,:)-xc(2,:))/norm(xc(1,:)-xc(3,:));
    sinth_0 =  norm(cross(([xc(2,:) 0]-[xc(1,:) 0]),([xc(3,:) 0]-[xc(1,:) 0])))/...
        norm(xc(1,:)-xc(2,:))/norm(xc(1,:)-xc(3,:));
    th_0 =  atan2(sinth_0,costh_0);
    rotIndx(p,5) = th_0;
end


end


