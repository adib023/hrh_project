function dp_dx =  cal_dp_dx(xCoord,indx,rotIndx, bndry, pbc)

nDof = 2*size(xCoord,1);
dp_dx = zeros(nDof,1);

for elem = 1:size(indx,1)
    n1 = indx(elem,1);
    n2 = indx(elem,2);
    k = indx(elem,3);
    xc = xCoord([n1,n2],:);
    l_0 = indx(elem,4);
    
    f_a =  getAxialForce(xc,k,l_0);
    
    if any(pbc(:,4) == n1)
        if any(bndry(:,1) == n2)
            
            id =  find(pbc(:,4) == n1);
            n1 = bndry(id,4);
        end
    end
    
    
    dp_dx([2*n1-1,2*n1,2*n2-1,2*n2],:) = dp_dx([2*n1-1,2*n1,2*n2-1,2*n2],:)+...
        f_a;
end
% 
dp_dx (pbc(:)) = 0;


for elem =  1:size(rotIndx,1)
    n1 = rotIndx(elem,1);
    n2 = rotIndx(elem,2);
    n3 = rotIndx(elem,3);
    xc = xCoord([n1 n2 n3],:);
    k =  rotIndx(elem,4);
    th_0 = rotIndx(elem,5);
    
    
    f_t =  getTorForce(xc,k,th_0);
%     
   [n1,n2,n3] = getCorPBCID(n1,n2,n3, bndry, pbc);
    
    dp_dx([2*n1-1,2*n1,2*n2-1,2*n2,2*n3-1,2*n3],:) = ...
        dp_dx([2*n1-1,2*n1,2*n2-1,2*n2,2*n3-1,2*n3],:)+ f_t;
    
end


end


function[n1,n2,n3] =  getCorPBCID(n1,n2,n3, bndry, pbc)


if any(bndry(:,1) == n1)
    if any(pbc(:,4) == n2) || any(pbc(:,4) == n3)
        if any(pbc(:,4) == n2)
            id = find(pbc(:,4) == n2);
            n2 = bndry(id,4);
        else
            id = find(pbc(:,4) == n3);
            n3 = bndry(id,4);
        end
    end
end

if any(bndry(:,4) == n1)
    if any(pbc(:,1) == n2) || any(pbc(:,1) == n3)
        if any(pbc(:,1) == n2)
            id = find(pbc(:,1) == n2);
            n2 = bndry(id,1);
        else
            id = find(pbc(:,1) == n3);
            n3 = bndry(id,1);
        end
    end
end


end