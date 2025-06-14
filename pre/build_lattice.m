function build_lattice(theta, pattern,N1, N2_list, k_list, kt_list, zeta, TB_disp)

N2 = sum(N2_list);
nNode_all = 0;
start = [0 0];
xCoord = [];
indx = [];
rotIndx = [];
ucID = [];


for n = 1:length(pattern)
    
    [a1,a2,start,nNode_all, xCoord1, indx1, rotIndx1,ucID1] =...
        lattice(theta,N1,N2_list(n),k_list(n),kt_list(n),start,pattern(n),nNode_all);
    size(rotIndx1,1)
    xCoord =  [xCoord;xCoord1];
    indx = [indx;indx1];
    rotIndx = [rotIndx;rotIndx1];
    ucID = [ucID ucID1];
end
size(rotIndx,1)
[indx, rotIndx] = ...
    interface(k_list(2),kt_list(2),pattern,N1, N2_list, xCoord, indx, rotIndx, ucID);
size(rotIndx,1)
rotIndx = fix_set_angle(xCoord,rotIndx);
[axialPBC, rotPBC] = define_PBC(N1,N2,xCoord,indx,rotIndx,ucID);
[disp_bndry, frcID] =  bndry_and_frc(N1,N2,ucID, TB_disp);
dmpr = damping(indx,axialPBC,zeta);


file =  strcat("../file/", pattern, '.txt' );
fig = pattern;

fid =  fopen(file, 'w');

numRot =  size(rotIndx,1);

write(file,a1);
write(file,a2);
write(file,[N1 N2_list]);
write(file,size(xCoord,1));
write(file,xCoord);
write(file,size(indx,1));
write(file,indx);
write(file,size(rotIndx,1));
write(file,rotIndx);
write(file,axialPBC);
write(file,rotPBC);
write(file,disp_bndry);
write(file,frcID);
write(file, dmpr)
write(file,ucID(:,:,1));
write(file,ucID(:,:,2));
write(file,ucID(:,:,3));
write(file,ucID(:,:,4));

fclose(fid);

plot_lattice(xCoord,indx,rotIndx,fig);


end