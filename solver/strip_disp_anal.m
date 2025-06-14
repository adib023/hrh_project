%%%%routine to determine the dispersion of strip

currentFolder =  pwd;

if currentFolder ~= "../solver"
    cd ../solver
end

clc
clear all
%close all

addpath("../pre")

%problemType = 1; %%% dispersion analysis for both
%problemType = 2; %%% modeshape plot for undeformed
 % problemType = 3; %%% modeshape plot for deformed
problemType = 4; %%% frequency Response analysis for both


fileName =  "../file/hrh.txt";
lattice =  read_lattice(fileName);

%%% determine a strips node

N1 = lattice.N1;
N2 = sum(lattice.N2_list);

id1 = lattice.id1;
id2 = lattice.id2;
id3 = lattice.id3;
id4 = lattice.id4;

stripNode = [];
rightNode = [];
leftNode = [];




for p = 1:N2
    stripNode = [stripNode;id1(N1-1,p);...
        id2(N1-1,p);id3(N1-1,p);id4(N1-1,p)];

    rightNode = [rightNode;id1(N1,p);...
        id2(N1,p);id3(N1,p);id4(N1,p)];

    leftNode =  [leftNode;id1(N1-2,p);...
        id2(N1-2,p);id3(N1-2,p);id4(N1-2,p)];
end


% stripNode(33)
% stripNode(36)

fixedNode = [];%lattice.presDispID;

strip_fix_check = ismember(stripNode, fixedNode);
strip_fix_id =  find(strip_fix_check);
strip_fix =  stripNode(strip_fix_id);



[W,Wint, left_rows, right_rows] = ...
    transformation_matrix_strip(lattice.nNode, ...
    strip_fix,stripNode, rightNode, leftNode);




if problemType  == 1

    xCoord =  lattice.xCoord;
    [K,f] = getStiffness_matrix(xCoord, lattice.indx, lattice.rotIndx,1,0);

    mu = -0.25*pi:0.001:0.25*pi;
    count = 0;

    freq = [];

    for m = mu
        count =  count + 1;
        Kred = ...
            getKred_strip(K,W,Wint,m, right_rows, left_rows);

        w2 = eig(Kred);
        w = sqrt(real(w2));
        w = sort(w);

        freq =  [freq w];
    end

    figure
    for p = 1: size(freq,1)
        plot(mu./pi,real(freq(p,:)), 'k-')
        hold on
    end


    xCoord =  read('../result/NRresult/def_xCoord.txt');
    [K,f] = getStiffness_matrix(xCoord, lattice.indx, lattice.rotIndx,1,0);


    mu = -0.25*pi:0.001:0.25*pi;
    count = 0;

    freq = [];

    for m = mu
        count =  count + 1;
        Kred = ...
            getKred_strip(K,W,Wint,m, right_rows, left_rows);

        w2 = eig(Kred);
        w = sqrt(real(w2));
        w = sort(w);

        freq =  [freq w];
    end

    for p = 1: size(freq,1)
        plot(mu./pi,real(freq(p,:)), 'r--')
        hold on
    end

    xlim([-0.25 0.25]);
    ylim([0.25 0.3])
    xlabel("$\mu/\pi$",'Interpreter','latex')
    ylabel("$\omega$",'Interpreter','latex')
    set(gca,'fontsize', 20)
    set(gca, 'fontname', "Times New Roman")


elseif problemType == 2
mu = 0;
     freqID = 42;
    xCoord =  lattice.xCoord;
    [K,f] = getStiffness_matrix(xCoord, lattice.indx, lattice.rotIndx,1,0);

    Kred = ...
        getKred_strip(K,W,Wint,mu, right_rows, left_rows);

    [A, w2] = eig(Kred,'vector');

    w = sqrt(real(w2));
    [w, ind] = sort(w);

    A = A(:,ind);
    % T = eye(length(stripNode));
    % T(:,strip_fix_id) = [];
    % A = T*A;

    plot_strip_mode(xCoord, lattice.indx, A, w, stripNode, problemType, freqID)
elseif problemType == 3

    mu = 0;
    freqID = 30;
    xCoord =  read('../result/NRresult/def_xCoord.txt');
    [K,f] = getStiffness_matrix(xCoord, lattice.indx, lattice.rotIndx,1,0);

    Kred = ...
        getKred_strip(K,W,Wint,mu, right_rows, left_rows);

    [A, w2] = eig(Kred,'vector');

    w = sqrt(real(w2));
    [w, ind] = sort(w);

    A = A(:,ind);
    % T = eye(length(stripNode));
    % T(:,strip_fix_id) = [];
    % A = T*A;

    plot_strip_mode(xCoord, lattice.indx, A, w, stripNode,problemType,freqID)

elseif problemType == 4

    %xCoord =  lattice.xCoord;
    xCoord =  read('../result/NRresult/def_xCoord.txt');
    [K,f] = getStiffness_matrix(xCoord, lattice.indx, lattice.rotIndx,1,0);

    Kred = ...
        getKred_strip(K,W,Wint,0, right_rows, left_rows);
    Mred =  eye(size(Kred));

    Fred =  zeros(size(Mred,1),1);

    Fred(34*2-1) = 1;
    Fred(35*2-1) = -1;
    
    U1_list = [];
    U2_list = [];
    for w =  0.26:0.000001:0.29
        U =  (-w^2*Mred + Kred)\Fred;
        U1_list = [U1_list norm(U(34*2-1:34*2))];
        U2_list = [U2_list norm(U(19*2-1:19*2))];
    end

    %figure
    plot(0.26:0.000001:0.29, 10*log10((U1_list./sqrt(2)).^2),'k--')
    hold on
    plot(0.26:0.000001:0.29, 10*log10((U2_list./sqrt(2)).^2),'r--')
    hold on
end




