function K_r = getRotK(xc,k_rot,th_0)

x1 = xc(1,1); y1 = xc(1,2);
x2 = xc(2,1); y2 = xc(2,2);
x3 = xc(3,1); y3 = xc(3,2);

costh =  dot((xc(2,:)-xc(1,:)),(xc(3,:)-xc(1,:)))/...
    norm(xc(1,:)-xc(2,:))/norm(xc(1,:)-xc(3,:));
sinth =  norm(cross(([xc(2,:) 0]-[xc(1,:) 0]),([xc(3,:) 0]-[xc(1,:) 0])))/...
    norm(xc(1,:)-xc(2,:))/norm(xc(1,:)-xc(3,:));

th =  atan2(sinth,costh);

l_12 = norm(xc(1,:)-xc(2,:));
l_13 = norm(xc(1,:)-xc(3,:));

Dx = zeros(3,2);
%%%%%% first derivative %%%%%%%
Dx(2,1) =  (y2-y1)/l_12^2;
Dx(2,2) =  -(x2-x1)/l_12^2;
Dx(3,1) =  -(y3-y1)/l_13^2;
Dx(3,2) =  (x3-x1)/l_13^2;
Dx(1,1) =  -Dx(2,1)-Dx(3,1);
Dx(1,2) =  -Dx(2,2)-Dx(3,2);

% Dx_debug = Dx * k_rot
% dx = debug_rotK(xc,k_rot,th_0)

DxDx = zeros(6);


DxDx(3,3) = -2*(x2-x1)*(y2-y1)/l_12^4;
DxDx(3,4) = ((x2-x1)^2-(y2-y1)^2)/l_12^4;
DxDx(3,1) = -DxDx(3,3);
DxDx(3,2) = -DxDx(3,4);


DxDx(4,3) = DxDx(3,4);
DxDx(4,4) = -DxDx(3,3);
DxDx(4,1) = DxDx(3,2);
DxDx(4,2) = -DxDx(4,4);


DxDx(5,5) = 2*(x3-x1)*(y3-y1)/l_13^4;
DxDx(5,6) = -((x3-x1)^2-(y3-y1)^2)/l_13^4;
DxDx(5,1) = -DxDx(5,5);
DxDx(5,2) = -DxDx(5,6);


DxDx(6,5) = DxDx(5,6);
DxDx(6,6) = -DxDx(5,5);
DxDx(6,1) = DxDx(5,2);
DxDx(6,2) = -DxDx(6,6);

DxDx(1,1) = (DxDx(3,3)+DxDx(5,5));
DxDx(1,2) = (DxDx(3,4)+DxDx(5,6));
DxDx(1,3) = -DxDx(3,3);
DxDx(1,4) = -DxDx(3,4);
DxDx(1,5) = -DxDx(5,5);
DxDx(1,6) = -DxDx(5,6);

DxDx(2,1) = DxDx(1,2);
DxDx(2,2) = DxDx(4,4)+DxDx(6,6);
DxDx(2,3) = DxDx(3,2);
DxDx(2,4) = DxDx(4,2);
DxDx(2,5) = DxDx(5,2);
DxDx(2,6) = DxDx(6,2);


DxDx = (th-th_0)*DxDx;

%%%%%%% derive the 6x6 stiffness matrix in the form of 
%%%%%%% \partial^2 P /\partial x_p \partial x_q = DxDx 
% 
DxDx(1,1) = DxDx(1,1)+Dx(1,1)*Dx(1,1);
DxDx(1,2) = DxDx(1,2)+Dx(1,1)*Dx(1,2);
DxDx(1,3) = DxDx(1,3)+Dx(1,1)*Dx(2,1);
DxDx(1,4) = DxDx(1,4)+Dx(1,1)*Dx(2,2);
DxDx(1,5) = DxDx(1,5)+Dx(1,1)*Dx(3,1);
DxDx(1,6) = DxDx(1,6)+Dx(1,1)*Dx(3,2);

DxDx(2,1) = DxDx(2,1)+Dx(1,2)*Dx(1,1);
DxDx(2,2) = DxDx(2,2)+Dx(1,2)*Dx(1,2);
DxDx(2,3) = DxDx(2,3)+Dx(1,2)*Dx(2,1);
DxDx(2,4) = DxDx(2,4)+Dx(1,2)*Dx(2,2);
DxDx(2,5) = DxDx(2,5)+Dx(1,2)*Dx(3,1);
DxDx(2,6) = DxDx(2,6)+Dx(1,2)*Dx(3,2);

DxDx(3,1) = DxDx(3,1)+Dx(2,1)*Dx(1,1);
DxDx(3,2) = DxDx(3,2)+Dx(2,1)*Dx(1,2);
DxDx(3,3) = DxDx(3,3)+Dx(2,1)*Dx(2,1);
DxDx(3,4) = DxDx(3,4)+Dx(2,1)*Dx(2,2);
DxDx(3,5) = DxDx(3,5)+Dx(2,1)*Dx(3,1);
DxDx(3,6) = DxDx(3,6)+Dx(2,1)*Dx(3,2);

DxDx(4,1) = DxDx(4,1)+Dx(2,2)*Dx(1,1);
DxDx(4,2) = DxDx(4,2)+Dx(2,2)*Dx(1,2);
DxDx(4,3) = DxDx(4,3)+Dx(2,2)*Dx(2,1);
DxDx(4,4) = DxDx(4,4)+Dx(2,2)*Dx(2,2);
DxDx(4,5) = DxDx(4,5)+Dx(2,2)*Dx(3,1);
DxDx(4,6) = DxDx(4,6)+Dx(2,2)*Dx(3,2);

DxDx(5,1) = DxDx(5,1)+Dx(3,1)*Dx(1,1);
DxDx(5,2) = DxDx(5,2)+Dx(3,1)*Dx(1,2);
DxDx(5,3) = DxDx(5,3)+Dx(3,1)*Dx(2,1);
DxDx(5,4) = DxDx(5,4)+Dx(3,1)*Dx(2,2);
DxDx(5,5) = DxDx(5,5)+Dx(3,1)*Dx(3,1);
DxDx(5,6) = DxDx(5,6)+Dx(3,1)*Dx(3,2);

DxDx(6,1) = DxDx(6,1)+Dx(3,2)*Dx(1,1);
DxDx(6,2) = DxDx(6,2)+Dx(3,2)*Dx(1,2);
DxDx(6,3) = DxDx(6,3)+Dx(3,2)*Dx(2,1);
DxDx(6,4) = DxDx(6,4)+Dx(3,2)*Dx(2,2);
DxDx(6,5) = DxDx(6,5)+Dx(3,2)*Dx(3,1);
DxDx(6,6) = DxDx(6,6)+Dx(3,2)*Dx(3,2);


% k_rot
% DxDx
K_r = k_rot*DxDx;

end

