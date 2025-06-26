function freq_response(M_global, K_global,forceNode, f_range, resNode, ...
    Wmat, freqFile, resultFile)

global DIM

fid1 = fopen(strcat('../Results/',freqFile),'w')';
fid1 = write(fid1, f_range);
fclose(fid1);

w_range = f_range * 2* pi;

extP = zeros(size(Wmat,1),1);

for p = 1:size(forceNode,1)
    for q = 1:DIM
        id =  (forceNode(p,1) - 1) *DIM + q;
        factr =  forceNode(p,q+1);
        extP(id) = factr;
    end
end

extP = Wmat' * extP;

fid2 =  fopen(strcat('../Results/',resultFile),'w');

for w =  w_range
    u = (-w^2*M_global + K_global)\extP;
    u = Wmat * u;
    u_norm_list =  [];
    for node = resNode
        u_norm = norm([u(3*node-2) u(3*node-1) u(3*node)]);
        u_norm_list = [u_norm_list u_norm/norm(extP)];
        
    end
    fid2 = write(fid2, u_norm_list);
end

fclose(fid2);

end
