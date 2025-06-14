currentFolder =  pwd;

if currentFolder ~= "../solver"
    cd ../solver
end

clc
clear all
%close all

addpath("../pre")

problemType = 1;
% problemType = 2;

fileName =  "../file/hrh.txt";
lattice =  read_lattice(fileName);

if problemType == 1
    xCoord =  lattice.xCoord;
else
    xCoord =  read('../result/NRresult/def_xCoord.txt');
end

[K,f] = getStiffness_matrix(xCoord, lattice.indx, lattice.rotIndx,1,0);

[K,~] = apply_PBC(xCoord,lattice.axialPBC,lattice.rotPBC,K,f,1,0);


M =  eye(size(K));

F =  zeros(size(M,1),1);

F(135*2-1) = 1;
F(134*2-1) = -1;


K([2*lattice.presDispID-1,2*lattice.presDispID],:)= [];
K(:,[2*lattice.presDispID-1,2*lattice.presDispID])= [];

M([2*lattice.presDispID-1,2*lattice.presDispID],:)= [];
M(:,[2*lattice.presDispID-1,2*lattice.presDispID])= [];

F([2*lattice.presDispID-1,2*lattice.presDispID])= [];

U1_list = [];
U2_list = [];
wList  = 0.26:0.00001:0.29;

for w =  wList
    U =  (-w^2*M + K)\F;

    T = eye(lattice.nNode*2);
    T(:,[2*lattice.presDispID-1,2*lattice.presDispID]) = [];
    U = T*U;

    U1_list = [U1_list norm(U(134*2-1:134*2))];
    U2_list = [U2_list norm(U(60*2-1:60*2))];
end


if problemType == 1
    figure
    plot(wList, 10*log10((U1_list./sqrt(2)).^2),'b-', 'LineWidth',1)
    hold on
    plot(wList, 10*log10((U2_list./sqrt(2)).^2),'r-', 'LineWidth',1)
    hold on

    set(gca,'fontname','times new roman')
set(gca,'fontsize',20)
xlabel('$\omega$','interpreter','latex')
ylabel('$|u/F|^2$(dB)','interpreter','latex')
else
    plot(wList, 10*log10((U1_list./sqrt(2)).^2),'b--', 'LineWidth',1)
    hold on
    plot(wList, 10*log10((U2_list./sqrt(2)).^2),'r--', 'LineWidth',1)
    hold on
end