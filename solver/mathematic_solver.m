close all
clear all
clc


ka = 1;
kt = 0.01;

dummy = 0;

th = -15* pi/ 180; 

a  = [2 * cos(th), 0 ];

x1 =  [0 0];
x2  = x1 + cos(th)/ cos(pi/6)*[cos( 5* pi/6) sin(5* pi/6)];
x3 = x2 + [0 1];
x4 = x3 + [cos(th) sin(th)];
x5 = x4 + [0 1];
x6 = x5 + [cos(pi - th) sin(pi - th)];
x7 = x6 + [0 1];
x8 = x7 + cos(th)/cos(pi/6) * [cos(pi/6) sin(pi/6)];

numNodeUC  = 8;
numUC =  4;

xCoord1 = [x1;x2;x3;x4;x5;x6;x7;x8];
xCoord = [];
nodeIDList = [];

for p = 1: numUC
    xCoord =  [xCoord; xCoord1 + (p-1)*a];
    nodeIDList =  [nodeIDList; (1:8) + (p-1)*numNodeUC;];
end

axialIndx = [];

for p = 1 : numUC
    for q = 1 : numNodeUC -1
        axialIndx = [axialIndx; [nodeIDList(p,q) nodeIDList(p,q+1) ka]];
    end
end

size(axialIndx)

for p = 1 : numUC-1
    
    axialIndx =  [axialIndx; [nodeIDList(p,1) nodeIDList(p+1,2) ka]];
    axialIndx =  [axialIndx; [nodeIDList(p,4) nodeIDList(p+1,3) ka]];
    axialIndx =  [axialIndx; [nodeIDList(p,5) nodeIDList(p+1,6) ka]];
    axialIndx =  [axialIndx; [nodeIDList(p,8) nodeIDList(p+1,7) ka]];

end

axialIndx(:,4) = dummy;
axialIndx(:,5) = dummy;


size(axialIndx)

rotIndx =  [nodeIDList(1,1) nodeIDList(2,2) nodeIDList(1,2);
            nodeIDList(1,2) nodeIDList(1,1) nodeIDList(1,3);
            nodeIDList(1,3) nodeIDList(1,2) nodeIDList(1,4);
            nodeIDList(1,4) nodeIDList(1,3) nodeIDList(2,3);
            nodeIDList(1,4) nodeIDList(1,5) nodeIDList(1,3);
            nodeIDList(1,4) nodeIDList(2,3) nodeIDList(1,5);
            nodeIDList(1,5) nodeIDList(1,6) nodeIDList(1,4);
            nodeIDList(1,5) nodeIDList(2,6) nodeIDList(1,6);
            nodeIDList(1,5) nodeIDList(1,4) nodeIDList(2,6);
            nodeIDList(1,6) nodeIDList(1,5) nodeIDList(1,7);
            nodeIDList(1,7) nodeIDList(1,6) nodeIDList(1,8);
            nodeIDList(1,8) nodeIDList(1,7) nodeIDList(2,7)];



for p = 2:numUC - 1
    
    rotIndx_p =  [nodeIDList(p,1) nodeIDList(p+1,2) nodeIDList(p,2);

            nodeIDList(p,2) nodeIDList(p,1) nodeIDList(p,3);
            nodeIDList(p,2) nodeIDList(p,3) nodeIDList(p-1,1);
            nodeIDList(p,2) nodeIDList(p-1,1) nodeIDList(p,1);

            nodeIDList(p,3) nodeIDList(p,2) nodeIDList(p,4);
            nodeIDList(p,3) nodeIDList(p-1,4) nodeIDList(p,2);
            nodeIDList(p,3) nodeIDList(p-1,4) nodeIDList(p,4);

            nodeIDList(p,4) nodeIDList(p+1,3) nodeIDList(p,3);
            nodeIDList(p,4) nodeIDList(p,5) nodeIDList(p,3);
            nodeIDList(p,4) nodeIDList(p+1,3) nodeIDList(p,5);

            nodeIDList(p,5) nodeIDList(p,6) nodeIDList(p,4);
            nodeIDList(p,5) nodeIDList(p,6) nodeIDList(p+1,6);
            nodeIDList(p,5) nodeIDList(p,4) nodeIDList(p+1,6);
            
            nodeIDList(p,6) nodeIDList(p,5) nodeIDList(p,7);
            nodeIDList(p,6) nodeIDList(p,7) nodeIDList(p-1,5);
            nodeIDList(p,6) nodeIDList(p,5) nodeIDList(p-1,5);
    
            nodeIDList(p,7) nodeIDList(p,6) nodeIDList(p,8);
            nodeIDList(p,7) nodeIDList(p-1,8) nodeIDList(p,6);
            nodeIDList(p,7) nodeIDList(p,8) nodeIDList(p-1,8);

            nodeIDList(p,8) nodeIDList(p,7) nodeIDList(p+1,7)];

    rotIndx = [rotIndx;rotIndx_p];

end

for p = 4
rotIndx =  [rotIndx;
            nodeIDList(p,2) nodeIDList(p,1) nodeIDList(p,3);
            nodeIDList(p,2) nodeIDList(p,3) nodeIDList(p-1,1);
            nodeIDList(p,2) nodeIDList(p-1,1) nodeIDList(p,1);

            nodeIDList(p,3) nodeIDList(p,2) nodeIDList(p,4);
            nodeIDList(p,3) nodeIDList(p-1,4) nodeIDList(p,2);
            nodeIDList(p,3) nodeIDList(p-1,4) nodeIDList(p,4);

            
            nodeIDList(p,4) nodeIDList(p,5) nodeIDList(p,3);
            
            nodeIDList(p,5) nodeIDList(p,6) nodeIDList(p,4);

            nodeIDList(p,6) nodeIDList(p,5) nodeIDList(p,7);
            nodeIDList(p,6) nodeIDList(p,7) nodeIDList(p-1,5);
            nodeIDList(p,6) nodeIDList(p,5) nodeIDList(p-1,5);
    
            nodeIDList(p,7) nodeIDList(p,6) nodeIDList(p,8);
            nodeIDList(p,7) nodeIDList(p-1,8) nodeIDList(p,6);
            nodeIDList(p,7) nodeIDList(p,8) nodeIDList(p-1,8)
    ];
end

rotIndx(:,4) =  kt;
rotIndx(:,5) =  dummy;
rotIndx(:,6) =  dummy;


rotIndx = fix_set_angle(xCoord,rotIndx)

%%%% define axial pbc now
noShift =  [0 0];
shift= 4 * a;  

axialPBC =  [25 2 ka noShift shift dummy;
    28 3 ka noShift shift dummy;
    29 6 ka noShift shift dummy;
    32 7 ka noShift shift dummy];


rotPBC =  [25 2 26 kt noShift shift noShift dummy;
           28 27 3 kt noShift noShift shift dummy;
           28 3 29 kt noShift shift noShift dummy;
           29 28 6 kt noShift noShift shift dummy;
           29 6 30 kt noShift shift noShift dummy;
           32 31 7 kt noShift noShift shift dummy;
           2 25 1 kt shift noShift shift dummy;
           2 3 25 kt shift shift noShift dummy;
           3 28 2 kt shift noShift shift dummy;
           3 4 28 kt shift shift noShift dummy;
           6 29 5 kt shift noShift shift dummy;
           6 7 29 kt shift shift noShift dummy;
           7 8 32 kt shift shift noShift dummy;
           7 32 6 kt shift noShift shift dummy
           ];







M =  eye(2*size(xCoord,1));
[K,f] = getStiffness_matrix(xCoord, axialIndx, rotIndx,1,0);
[K,~] = apply_PBC(xCoord,axialPBC,rotPBC,K,f,1,0);


fid =  fopen('stiffness_matrix.txt','w');

write('stiffness_matrix.txt', K)

fclose(fid);



[A, w2] =  eig(K,"vector");
[w2, ind] = sort(w2);
A = A(:,ind);

w =  real(sqrt(w2));

fid =  fopen('localized_mode.txt','w');

write('localized_mode.txt', w(34));
write('localized_mode.txt', A(:,34));

fclose(fid);


% figure
% plot(w,'o')
% 
% print('-dpng', '-r300', '../result/theory_result/frequency');
% 
% alpha = 1;
% 
% n = size(xCoord,1);
% for p = 1:length(w)
%     defCoord =  xCoord;
%     defCoord(:,1) =  defCoord(:,1) + alpha * A(1:2:2*n, p);
%     defCoord(:,2) =  defCoord(:,2) + alpha * A(2:2:2*n, p);
% 
%     figure
%     plot_lattice(defCoord, axialIndx, 'k-');
%     plot_lattice(xCoord, axialIndx, 'r--');
%     title(num2str(w(p)))
%     print('-dpng', '-r300', strcat('../result/theory_result/modeshape/',num2str(p)))
% end

