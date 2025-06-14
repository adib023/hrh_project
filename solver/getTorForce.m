function f_t = getTorForce(xc,k_rot,th_0)

f_t = zeros(6,1);

x1 =  xc(1,1); y1 =  xc(1,2);
x2 =  xc(2,1); y2 =  xc(2,2);
x3 =  xc(3,1); y3 =  xc(3,2);

costh =  dot((xc(2,:)-xc(1,:)),(xc(3,:)-xc(1,:)))/...
    norm(xc(1,:)-xc(2,:))/norm(xc(1,:)-xc(3,:));
sinth =  norm(cross(([xc(2,:) 0]-[xc(1,:) 0]),([xc(3,:) 0]-[xc(1,:) 0])))/...
    norm(xc(1,:)-xc(2,:))/norm(xc(1,:)-xc(3,:));

th =  atan2(sinth,costh);

l12 = norm(xc(1,:)-xc(2,:));
l13 = norm(xc(1,:)-xc(3,:));


f_t(3) = (y2-y1)/l12/l12;
f_t(4) = -(x2-x1)/l12/l12;
f_t(5) = -(y3-y1)/l13/l13;
f_t(6) = (x3-x1)/l13/l13;
f_t(1) = -f_t(3)-f_t(5);
f_t(2) = -f_t(4)-f_t(6);

f_t = f_t*(th-th_0) * k_rot;
end