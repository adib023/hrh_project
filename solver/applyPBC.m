function [K] = applyPBC(xCoord,bndryID,pbcID,K,mu)



for p = 1: size(bndryID,1)
    nb = bndryID(p,1);
    n =  pbcID(p,1);
    a = xCoord(n,1:2) - xCoord(nb,1:2);
    
    K(:,2*n-1:2*n) = K(:,2*n-1:2*n)*exp(i*dot(mu,a));
    K(:,2*nb-1:2*nb)= K(:,2*nb-1:2*nb)+ K(:,2*n-1:2*n);
    
    nb = bndryID(p,2);
    n =  pbcID(p,2);
    a = xCoord(n,1:2) - xCoord(nb,1:2);
    
    K(:,2*n-1:2*n) = K(:,2*n-1:2*n)*exp(i*dot(mu,a));
    K(:,2*nb-1:2*nb)= K(:,2*nb-1:2*nb)+ K(:,2*n-1:2*n);
    
    nb = bndryID(p,3);
    n =  pbcID(p,3);
    a = xCoord(n,1:2) - xCoord(nb,1:2);
    
    K(:,2*n-1:2*n) = K(:,2*n-1:2*n)*exp(i*dot(mu,a));
    K(:,2*nb-1:2*nb)= K(:,2*nb-1:2*nb)+ K(:,2*n-1:2*n);
    
    nb = bndryID(p,4);
    n =  pbcID(p,4);
    a = xCoord(n,1:2) - xCoord(nb,1:2);
    
    K(:,2*n-1:2*n) = K(:,2*n-1:2*n)*exp(i*dot(mu,a));
    K(:,2*nb-1:2*nb)= K(:,2*nb-1:2*nb)+ K(:,2*n-1:2*n);
end


end