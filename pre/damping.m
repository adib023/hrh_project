function dmpr = damping(indx,axialPBC,zeta)

dmpr1(:,1) = indx(:,1);
dmpr1(:,2) = indx(:,2);

for n = 1:size(indx,1)
    k = indx(n,3);
    dmpr1(n,3) = zeta*2*sqrt(k);
end

dmpr2(:,1) =  axialPBC(:,1);
dmpr2(:,2) =  axialPBC(:,2);

for n = 1:size(axialPBC,1)
    k = axialPBC(n,3);
    dmpr2(n,3) = zeta*2*sqrt(k);
end

dmpr = [dmpr1;dmpr2];

end