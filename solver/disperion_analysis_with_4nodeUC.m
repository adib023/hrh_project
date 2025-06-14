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
pattern = 'h';

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

kList = [OX;XM;MY;YO];


centerUCNOde =  ucID(2,2,:);
target  =  [2*centerUCNOde-1 2*centerUCNOde];
target = sort(target);
a1(1)
wList = [];
for p = 1:size(kList,1)
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

    w2 = eig(Kred);
    w2 = sort(w2);
    w  = sqrt(w2);
    w = real(w);
    wList = [wList w];
end


figure

if pattern ==  'r'
    line =  'r-';
else
    line = 'k-';
end
for n = 1:8
    plot(wList(n,:),line);
    hold on
end

xticks([0 100 200 300 400])
%xticklabels({'\Gamma','X','M','Y','\Gamma'})
xticklabels({'\Gamma','X','M','Y','\Gamma'})
ylabel("$\omega$","interpreter","latex","fontsize",20);
set(gca,'fontsize', 20);
set(gca,'fontname', 'times new roman')
ylim([0 1])
%xlim([300 400])
% yline(0.27073,'g--', LineWidth=2)



