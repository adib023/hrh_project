function [M,K] = M_K_mat(lattice,mu)

M = zeros(2*lattice.nNode);
K = zeros(2*lattice.nNode);

for n =  1 : lattice.nNode
    m =  1;
    M(2*n-1,2*n-1) = m;
    M(2*n,2*n) = m;
end

for elem =  1:lattice.nElem
    n1 = lattice.indx(elem,1);
    n2 = lattice.indx(elem,2);
    
    n1x =  2*n1- 1;
    n2x =  2*n2- 1;
    
    n1y =  2*n1;
    n2y =  2*n2;
    
    
    K(n1x, n1x) =  K(n1x, n1x) + lattice.indx(elem,3);
    K(n1x, n2x) =  K(n1x, n2x) - lattice.indx(elem,3);
    K(n2x, n1x) =  K(n2x, n1x) - lattice.indx(elem,3);
    K(n2x, n2x) =  K(n2x, n2x) + lattice.indx(elem,3);
    
    K(n1y, n1y) =  K(n1y, n1y) + lattice.indx(elem,3);
    K(n1y, n2y) =  K(n1y, n2y) - lattice.indx(elem,3);
    K(n2y, n1y) =  K(n2y, n1y) - lattice.indx(elem,3);
    K(n2y, n2y) =  K(n2y, n2y) + lattice.indx(elem,3);
    
end

for elem =  1: lattice.nPBC
    n1 =  lattice.pbc(elem,1);
    n2 =  lattice.pbc(elem,2);
    k = lattice.pbc(elem,3);
    
    n1x = 2*n1-1;
    n2x = 2*n2-1;
    
    
    ax =  lattice.ax;
    N1 = lattice.N1;
    
    
    K(n1x, n1x) =  K(n1x, n1x) + k;
    K(n1x, n2x) =  K(n1x, n2x) - k*exp(-i*N1*dot(ax,mu));
    K(n2x, n1x) =  K(n2x, n1x) - k*exp(-i*N1*dot(ax,mu));
    K(n2x, n2x) =  K(n2x, n2x) + k;
    
end 


end 