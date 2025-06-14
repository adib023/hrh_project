function [OA,AB,BC,CO] = getIBZ_re(a1,a2)

O = [0 0];

D =  [a1' a2'];
R =  D^-1; %%%%%reciprocal lattice vector matricx R = [r1;r2]

r1 = pi*R(1,:)'
r2 = pi*R(2,:)'

r12 =  r1+r2

theta_r1r2 =  acos(dot(r1,r2)/norm(r1)/norm(r2));
theta_r1r12 =  acos(dot(r1,r12)/norm(r1)/norm(r12));
theta_r2r12 =  acos(dot(r2,r12)/norm(r2)/norm(r12));

theta_r1r2_de = theta_r1r2 * 180/pi
theta_r1r12_de = theta_r1r12 * 180/pi
theta_r2r12_de = theta_r2r12 * 180/pi


A =  r1 + r12;
C =  r2;

B =  A + 2*(r12 - A);

figure
plot([O(1) A(1)],[O(2) A(2)],"b")
hold on
plot([A(1) B(1)],[A(2) B(2)],"r")
hold on
plot([B(1) C(1)],[B(2) C(2)],"k")
hold on
plot([C(1) O(1)],[C(2) O(2)],"g")
hold on


end