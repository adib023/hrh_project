function disp_anal_old()
[a1,a2,xCoord,indx,rotIndx] =  read("../file/hexaUC.txt");

D =  [a1' a2'];
B = D^-1;

b1 = B(1,:)';
b2 = B(2,:)';

OA =  pi*b1*linspace(0,1,100);
AB =  pi*b1+pi*b2*linspace(0.01,1,100);
BO =  pi*(b1+b2)*linspace(0.99,0,100);

k_vectors = [OA AB BO]

% ka = 1;
% kt = 0.1;
% 
% for n = 1: size(k_vectors,2)
%     w = hexa_disp_anal(xCoord,indx,rotIndx,ka,kt,k_vectors(:,n));
%     wList(:,n) = w;   
% end
% 
% 
% figure 
% for n = 1: size(k_vectors,2)
%     wSize = length(wList(:,n));
%     plot(n*ones(wSize,1),wList(:,n),"k.","markersize",5)
%     hold on
% end
% xlim([1 n]);



end