%%% this code is written on 5/8/2025
%%% this code is for doing dispersion analysis of different re-entrant
%%% with diffrent angle
%%% written by Adib Rahman - Texas A&M University

clc
clear all

addpath('../pre')



N1 =  3;
N2 =  3;
nNode_all = 0;
start = [0 0];

theta = -15;
pattern = 'r';

k = 1;
kt = 0.01;


[a1,a2,start,nNode_all, xCoord, indx, rotIndx,ucID] =...
    lattice(theta,N1,N2,k,kt,start,pattern,nNode_all);

rotIndx = fix_set_angle(xCoord,rotIndx);

%plot_lattice_pre(xCoord,indx,rotIndx,'disp');

[K,~] = getStiffness_matrix(xCoord, indx, rotIndx,1,0);

b1 = 1/a1(1)*pi;
b2 = 1/a2(2)*pi;

OX = [linspace(0,b1,100)'*a1(1) zeros(100,1)];
XM = [b1*a1(1)*ones(100,1) linspace(0,b2,100)'*a2(2)];
MY = [linspace(b1,0,100)'*a1(1) b2*a2(2)*ones(100,1)];
YO = [zeros(100,1) linspace(b2,0,100)'*a2(2)];

kList =[OX;XM;MY;YO];


centerUCNOde =  ucID(2,2,:);
target  =  [2*centerUCNOde-1 2*centerUCNOde];
target = sort(target);
a1(1)
wList = [];
for p =  390
    mux = kList(p,1);
    muy = kList(p,2);
    Kred = K(target,:);

    for n1 =  1:3
        for n2 = 1:3
            if n1 == 2 && n2 == 2
                continue
            else
                node1 =  ucID(n1,n2,1);
                node2 =  ucID(n1,n2,2);
                node3 =  ucID(n1,n2,3);
                node4 =  ucID(n1,n2,4);

                pn =  [2*[node1 node2 node3 node4]-1 2*[node1 node2 node3 node4]];
                pn = sort(pn);

                Kred(:,target) = Kred(:,target) + Kred(:,pn)*exp(i*dot([n1-2 n2-2],[mux,muy]));   


            end
        end
    end
    Kred =  Kred(:,target);
    [A,w2] = eig(Kred,"vector");
    [w2,ind] = sort(w2); 
    A = A(:,ind);
    w  = sqrt(w2);
    w = real(w);
    wList = [wList w];
end

% length(w)
% 
% for q = 1:length(w)
%     q
%         eigV = A(:,q);
%         dot_sum = 0;
%         cross_sum = 0;
%         for r = 1:4
%             nodeU  = eigV(2*r-1:2*r)';
%             dot_sum  = dot_sum + abs(dot([mux muy 0], [nodeU 0]))
%             cross_sum = cross_sum + norm(cross([mux muy 0], [nodeU 0]))
%         end
% 
%         % if dot_sum > cross_sum
%         %     wList_p = [wList_p;[p w(q)]]
%         % else
%         %     wList_s = [wList_s;[p w(q)]]
%         % end
% end


for r = 1:3
    eigV = A(:,r);
figure

for p = 1:size(indx,1)
    n1 = indx(p,1);
    n2 = indx(p,2);
    
    
    check = ismember([n1 n2], target);
    if sum(check) < 2
        continue
    end
    

    t1 =  find(target(:) == n1);
    t2 =  find(target(:) == n2);
    xc1 = xCoord(n1,:);
    xc2 = xCoord(n2,:);

    plot([xc1(1) xc2(1)], [xc1(2) xc2(2)], 'r--');
    hold on

    
    xc1 = xCoord(n1,:)+eigV(t1*2-1:t1*2)';
    xc2 = xCoord(n2,:)+eigV(t2*2-1:t2*2)';

    plot([xc1(1) xc2(1)], [xc1(2) xc2(2)], 'k-');
    hold on
end
end