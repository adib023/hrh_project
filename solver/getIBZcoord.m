function [OA,AB,BC,CD,DO]=getIBZcoord(a1,a2)
O = [0;0];

% D =  [a1' a2'];
% B = D^-1;
% b1 = pi*B(1,:)';
% b2 = pi*B(2,:)';


%%%%%% for re entrant having corner angle of 15 degree
O = [0;0];
A = [4.38817;0];
B = [0;3.36715];
C = [1.13581;4.42365];
D = [0;4.23863];

OA = O + (A-O)*linspace(0,1,1000);
AB = A + (B-A)*linspace(0.0001,1,1000);
BC = B + (C-B)*linspace(0.0001,1,1000);
CD = C + (D-C)*linspace(0.0001,1,1000);
DO = D + (O-D)*linspace(0.0001,1,1000);


%%%%% for reg hexagon
% B = b1;
% normOA =  norm(b1)/cos(pi/6);
% A = [normOA;0];
% 
% OA = O + (A-O)*linspace(0,1,100);
% AB = A + (B-A)*linspace(0,1,100);
% BO = B + (O-B)*linspace(0,1,100);



end