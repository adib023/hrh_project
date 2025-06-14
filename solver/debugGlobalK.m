clc
clear all
close all
addpath("../pre")

filename = "../file/re_hexa_re_pbc.txt";
nDofNode = 2;
lattice = read_lattice(filename);

tol =  1e-7;

debugK = zeros(lattice.nNode*2);

for p = 1:lattice.nNode*2
    if mod(p,2) == 0
        n = p/2;
        s_xCoord = lattice.xCoord;
        s_xCoord(n,2)= s_xCoord(n,2) + tol; 
    else
        n = (p+1)/2;
        s_xCoord = lattice.xCoord;
        s_xCoord(n,1)= s_xCoord(n,1) + tol;
    end
    
    for q = 1:lattice.nNode*2
        df = first_derivative(lattice.xCoord,tol,q,lattice.indx,lattice.rotIndx);
        df1= first_derivative(s_xCoord,tol,q,lattice.indx,lattice.rotIndx);
        
        debugK(p,q) = (df1-df)/tol;
        debugK(q,p) = (df1-df)/tol;
    end
    
end
% 
% debugK
K = getStiffness_matrix(lattice.xCoord,lattice.indx,lattice.rotIndx);

K-debugK



function df = first_derivative(xCoord,tol,dof,indx,rotIndx)
    if mod(dof,2) == 0
        n = dof/2;
        f_xCoord = xCoord;
        f_xCoord(n,2)= f_xCoord(n,2) + tol; 
    else
        n = (dof+1)/2;
        f_xCoord = xCoord;
        f_xCoord(n,1)= f_xCoord(n,1) + tol;
    end
    
    E = getTotalE(xCoord,indx,rotIndx);
    E1 = getTotalE(f_xCoord,indx,rotIndx);
    
    df = (E1-E)/tol;
end