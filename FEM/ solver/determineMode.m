function [f,V] = determineMode(M_global, K_global, f, numMode, Wmat, freqFile, modeFile)




[V, w2] =  eigs(K_global, M_global, numMode, (2*pi*f)^2 );

w2= diag(w2);

[w2, ind] =  sort(w2);
V = V(:,ind);



w = real(sqrt(w2));

f = w/2/pi;

V  = Wmat * V;

fid1 = fopen(strcat('../Results/',freqFile),'w')';
fid2 =  fopen(strcat('../Results/',modeFile),'w');

fid1 = write(fid1, f);
fid2 = write(fid2, V);

fclose(fid1);
fclose(fid2);



end
