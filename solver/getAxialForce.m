function f_a =  getAxialForce(xc,k,l_0)

f_a =  zeros(4,1);

xv1 = xc(1,:);
xv2 = xc(2,:);

dx =  xv1(1) - xv2(1);
dy =  xv1(2) - xv2(2);
l12 =  norm(xv1-xv2);

%f_a  = [dp_a/d_x1, dp_a/d_y1 dp_a/d_x2 dp_a/d_y2]

f_a(1) = k*dx/l12;
f_a(2) = k*dy/l12;
f_a(3) = -f_a(1);
f_a(4) = -f_a(2);

f_a = (l12-l_0)*f_a;

end