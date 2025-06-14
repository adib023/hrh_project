function dx = debug_rotK(xc,k,th_0)

% n1 = rotIndx(1,1);
% n2 = rotIndx(1,2);
% n3 = rotIndx(1,3);
% 
% xc(1,1:2) = xCoord(n1,1:2);
% xc(2,1:2) = xCoord(n2,1:2);
% xc(3,1:2) = xCoord(n3,1:2);
% 
% k = rotIndx(1,4);
% [n1 n2 n3]

vec12 = xc(2,1:2)-xc(1,1:2);
vec13 = xc(3,1:2)-xc(1,1:2);



% costh = dot(vec12,vec13)/norm(vec12)/norm(vec13);
% sinth = norm(cross([vec12 0],[vec13 0]))/norm(vec12)/norm(vec13);
% th_0 = atan2(sinth,costh);

[E,th] = cal_energy(xc,k,th_0);

tol =  1e-6;

for p =1:3
    for q = 1:2
        xc(p,q) = xc(p,q) + tol;
        [E1,th] = cal_energy(xc,k,th_0);
        xc(p,q) = xc(p,q) - tol;
        dx(p,q) = (E1-E)/tol/(th-th_0);
    end 
end



end


function [E,th] = cal_energy(xc,k,th_0)

vec12 = xc(2,1:2)-xc(1,1:2);
vec13 = xc(3,1:2)-xc(1,1:2);



costh = dot(vec12,vec13)/norm(vec12)/norm(vec13);
sinth = norm(cross([vec12 0],[vec13 0]))/norm(vec12)/norm(vec13);

th = atan2(sinth,costh);

E =  0.5*k*(th-th_0)^2;

end