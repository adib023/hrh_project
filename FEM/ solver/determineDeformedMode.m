function determineDeformedMode(M_global, K_global, f, numMode, Wmat)

[V, w2] =  eigs(inv(M_global)* K_global, numMode, (2*pi*f)^2 );

w2= diag(w2);
w = real(sqrt(w2));

V  = Wmat * V;

fid1 = fopen('../Results/defFreq.txt','w')';
fid2 =  fopen("../Results/defMode.txt",'w');

fid1 = write(fid1, w);
fid2 = write(fid2, V);

fclose(fid1);
fclose(fid2);



end
