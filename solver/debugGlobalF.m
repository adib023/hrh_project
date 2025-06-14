clc
clear all
close all
addpath("../pre")

filename = "../file/re_hexa_re_pbc.txt";
nDofNode = 2;
lattice = read_lattice(filename);
E = getTotalE(lattice.xCoord,lattice.indx,lattice.rotIndx)
tol =  1e-10;

debugF = zeros(lattice.nNode*2,1);

for p = 1:lattice.nNode*2
   if mod(p,2) == 0
        n = p/2;
        f_xCoord = lattice.xCoord;
        f_xCoord(n,2)= f_xCoord(n,2) + tol; 
    else
        n = (p+1)/2;
        f_xCoord = lattice.xCoord;
        f_xCoord(n,1)= f_xCoord(n,1) + tol;
    end
    
  
    E1 = getTotalE(f_xCoord,lattice.indx,lattice.rotIndx);
    
    debugF(p) = (E1-E)/tol;
end


debugF
F =  cal_dp_dx(lattice.xCoord,lattice.indx,lattice.rotIndx,lattice.bndry, lattice.pbc)


err =debugF-F