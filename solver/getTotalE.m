function E = getTotalE(xCoord,indx,rotIndx)

E = 0;

for elem = 1:size(indx,1)
    n = indx(elem,1:2);
    k = indx(elem,3);
    l_0 = indx(elem,4);
    xc =  xCoord(n,:);
    Ea = getAxialE(xc,k,l_0);
    E = E + Ea;
end

for elem = 1:size(rotIndx,1)
    n = rotIndx(elem,1:3);
    k = rotIndx(elem,4);
    th_0 = rotIndx(elem,5);
    xc =  xCoord(n,:);
    Et = getTorE(xc,k,th_0);
    E = E + Et;
end
end


function Et = getTorE(xc,k,th_0)

vec12 = xc(2,1:2)-xc(1,1:2);
vec13 = xc(3,1:2)-xc(1,1:2);
costh = dot(vec12,vec13)/norm(vec12)/norm(vec13);
sinth = norm(cross([vec12 0],[vec13 0]))/norm(vec12)/norm(vec13);

th = atan2(sinth,costh);

Et =  0.5*k*(th-th_0)^2;

end

function Ea = getAxialE(xc,k,l_0)


xv1 = xc(1,:);
xv2 = xc(2,:);
l =  norm(xv1-xv2);

Ea =  0.5*k*(l-l_0)^2;

end