function eigenSolver(lattice)

nDofNode = 2;

%%%%%% define mass matrix %%%%%%

% M = eye(nDofNode*lattice.nNode);

%%%%%%%%%% mu = [mux, muy]%%%%%%%%%%
[K,f] = getStiffness_matrix(lattice.xCoord, lattice.indx, lattice.rotIndx,1,0);
[K,~] = apply_PBC(lattice.xCoord,lattice.axialPBC,lattice.rotPBC,K,f,1,0);

fid = fopen('stiffness.txt','w');

write('stiffness.txt',K);

fclose(fid);

K([2*lattice.presDispID-1,2*lattice.presDispID],:)= [];
K(:,[2*lattice.presDispID-1,2*lattice.presDispID])= [];


[A,w2] = eig(K,"vector");
[w2,ind] = sort(w2);
A = A(:,ind);

T = eye(lattice.nNode*2);
T(:,[2*lattice.presDispID-1,2*lattice.presDispID]) = [];
A = T*A;


w = sqrt(w2);

length(w)
filename = "../result/eigenResult.txt";

f = fopen(filename, "w");

write(filename,w);
write(filename,real(A));
write(filename,imag(A));


fclose(f);

end 
