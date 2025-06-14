function u = setDispForPBC(bndryID,pbcID,u)

for p = 1: size(bndryID,1)
    
    nb = bndryID(p,1);
    n =  pbcID(p,3);
    u(n) = u(nb);
    
    nb = bndryID(p,2);
    n =  pbcID(p,4);
    u(n) = u(nb);
    
    nb = bndryID(p,3);
    n =  pbcID(p,1);
    u(n) = u(nb);
    
    nb = bndryID(p,4);
    n =  pbcID(p,2);
    u(n) = u(nb);
end

end