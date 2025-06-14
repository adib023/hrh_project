function K_a =  getAxialSpringK(xc,k,l_0)


xv1 = xc(1,:);
xv2 = xc(2,:);

dx = xv1(1) - xv2(1);
dy = xv1(2) - xv2(2);

l = norm(xv1-xv2);

dxdx =  zeros(4);

dxdx(1,1) = dy^2/l^3;
dxdx(1,2) = -dx*dy/l^3;
dxdx(1,3) = -dxdx(1,1);
dxdx(1,4) = dx*dy/l^3;

dxdx(2,1) = dxdx(1,2);
dxdx(2,2) = dx^2/l^3;
dxdx(2,3) = dx*dy/l^3;
dxdx(2,4) = -dxdx(2,2);

dxdx(3,1) = dxdx(1,3);
dxdx(3,2) = dxdx(2,3);
dxdx(3,3) = dxdx(1,1);
dxdx(3,4) = -dx*dy/l^3;

dxdx(4,1) = dxdx(1,4);
dxdx(4,2) = dxdx(2,4);
dxdx(4,3) = dxdx(3,4);
dxdx(4,4) = dxdx(2,2);

dxdx = (l-l_0)*k*dxdx;

K_a =  zeros(4);

K_a(1,1) = k*dx*dx/l/l;
K_a(1,2) = k*dx*dy/l/l;
K_a(1,3) = -K_a(1,1);
K_a(1,4) = -K_a(1,2);

K_a(2,1) = K_a(1,2);
K_a(2,2) = k*dy*dy/l/l;
K_a(2,3) = -K_a(1,2);
K_a(2,4) = -K_a(2,2);

K_a(3,1) = K_a(1,3);
K_a(3,2) = K_a(2,3);
K_a(3,3) = K_a(1,1);
K_a(3,4) = K_a(1,2);

K_a(4,1) = K_a(1,4);
K_a(4,2) = K_a(2,4);
K_a(4,3) = K_a(3,4);
K_a(4,4) = K_a(2,2);

K_a = K_a + dxdx;

end