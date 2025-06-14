function [errf,errK] = axialElemDebug(xc,u,k,d_0)

xc = xc + u;

f =  getAxialForce(xc,k,d_0);
K =  getAxialSpringK(xc,k,d_0);
TOL = 1e-6;
count = 0;

E = getE(xc,k,d_0);

for p = 1:size(xc,1)
    for q = 1:2
        count = count + 1;
        d_xc =  xc;
        d_xc(p,q) = xc(p,q)+ TOL;
        E1 = getE(d_xc,k,d_0);
        f1 = getAxialForce(d_xc,k,d_0);
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


xv1 = xc(1,:);
xv2 = xc(2,:);
d =  norm(xv1-xv2);

E =  0.5*k*(d-d_0)^2;

end