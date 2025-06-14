function [errf,errK] = torElemDebug(xc,u,k,d_0)
xc = xc + u;

f =  getTorForce(xc,k,d_0);
K =  getRotK(xc,k,d_0);
TOL = 1e-6;
count = 0;

E = getE(xc,k,d_0);

for p = 1:size(xc,1)
    for q = 1:2
        count = count + 1;
        d_xc =  xc;
        d_xc(p,q) = xc(p,q)+ TOL;
        E1 = getE(d_xc,k,d_0);
        f1 = getTorForce(d_xc,k,d_0);
        d_f(count,:) = (E1-E)/TOL;
        d_K(count,:) = (f1-f)/TOL;
    end
end
d_f
f
errf = d_f- f;
errK = d_K - K;
end


function E = getE(xc,k,d_0)

vec12 = xc(2,1:2)-xc(1,1:2);
vec13 = xc(3,1:2)-xc(1,1:2);

costh = dot(vec12,vec13)/norm(vec12)/norm(vec13);
sinth = norm(cross([vec12 0],[vec13 0]))/norm(vec12)/norm(vec13);

d = atan2(sinth,costh);

E =  0.5*k*(d-d_0)^2;

end