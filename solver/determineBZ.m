currentFolder = pwd;
if currentFolder ~= "../solver"
    cd ../solver
end

lattice = read_lattice('../file/h.txt');

FBZline = getFBZ(lattice);

[K,~] = getStiffness_matrix(lattice.xCoord, lattice.indx, lattice.rotIndx,1,0);

id1 = lattice.id1;
id2 = lattice.id2;
id3 = lattice.id3;
id4 = lattice.id4;


freqList = [];
for p  = 1:size(FBZline,1)
    mu =  FBZline(p,:);
    redK = getRedDispK(K,mu, id1,id2,id3,id4);
    w2 =  eig(redK);
    f = sort(real(sqrt((w2))));
    freqList = [freqList f];
end
size(freqList)
figure
for p = 1:size(freqList,1)
    plot(freqList(p,:),"k.")
    hold on
end





function redK = getRedDispK(K,mu, id1,id2,id3,id4)

redK = zeros(4*2);


n1 = id1(2,2);
n2 = id2(2,2);
n3 = id3(2,2);
n4 = id4(2,2);

K = K([n1*2-1 n1*2 n2*2-1 n2*2 n3*2-1 n3*2 n4*2-1 n4*2],:);

for p = 1:3
    for q = 1:3
        m =  [p q] - [2 2];
        n1 = id1(p,q);
        n2 = id2(p,q);
        n3 = id3(p,q);
        n4 = id4(p,q);

        redK = redK + K(:,[n1*2-1 n1*2 n2*2-1 n2*2 n3*2-1 n3*2 n4*2-1 n4*2]) * exp(i*dot(m,mu)); 
    end
end


end
