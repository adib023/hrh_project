%%%%% finite difference debug
clc
clear all
close all
TOL =  1e-6;

filename = "../file/re_hexa_re_pbc.txt";
lattice = read_lattice(filename);

nDof  =  lattice.nNode *2;
 
x = zeros(nDof,1);

% x(2*lattice.presDispID -1) = lattice.presDispX;
% x(2*lattice.presDispID) = lattice.presDispY;

def_xCoord =  lattice.xCoord;

for n =  1:size(def_xCoord,1)
    p = n;
%     if any(lattice.pbc(:) == n)
%         [row,col] = find(lattice.pbc ==  n);
%         p = lattice.bndry(row, col);
%     end
    
    def_xCoord(n,1) = def_xCoord(n,1) + x(2*p-1);
    def_xCoord(n,2) = def_xCoord(n,2) + x(2*p);
end


f = cal_dp_dx(def_xCoord,lattice.indx,lattice.rotIndx,lattice.bndry,lattice.pbc);   
K = getStiffness_matrix(lattice.xCoord, lattice.indx, lattice.rotIndx, def_xCoord);

%K = applyPBC(def_xCoord,lattice.bndry,lattice.pbc,K,[0 0]);

%elimID = [lattice.pbc(:)' lattice.presDispID'];

% K([2*elimID-1, 2*elimID],:) = [];
% K(:,[2*elimID-1, 2*elimID]) = [];
% f([2*elimID-1, 2*elimID],:) = [];

count = 0;
K_deb =  zeros(size(K));
for n =  1:size(lattice.xCoord,1)
%     if any(lattice.pbc(:) == n)
%         continue
%     elseif any(lattice.presDispID == n)
%         continue
%     end
    for m = 1:2
        count = count + 1;
        deb_xc = def_xCoord;
        deb_xc(n,m) = deb_xc(n,m)+TOL;
        f1 = cal_dp_dx(deb_xc,lattice.indx,lattice.rotIndx,lattice.bndry,lattice.pbc);
%         f1([2*elimID-1, 2*elimID],:) = [];
        df =  (f1-f)/TOL;
        K_deb(count,:) = df;  
    end
end




% for p = 1:size(K,1)
%     for q = 1:size(K,2)
%         err(p,q) = (K_deb(p,q) - K(p,q))/K(p,q)*100;
%     end
% end

K
K_deb
err = K- K_deb