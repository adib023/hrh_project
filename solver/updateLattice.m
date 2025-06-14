function [def_xCoord, rotIndx] =  updateLattice(x,xCoord,old_rotIndx)

%%%%% here xCoord means the initial xCoord before updating latest
%%%%% deformation

def_xCoord =  xCoord;

for n =  1:size(def_xCoord,1)
    def_xCoord(n,1) = def_xCoord(n,1) + x(2*n-1);
    def_xCoord(n,2) = def_xCoord(n,2) + x(2*n);
end

for elem = 1:size(old_rotIndx,1)
    rotIndx = old_rotIndx;
    n =  rotIndx(elem,1:3);
    xc =  xCoord(n,:);
    eta = cross(([xc(2,:) 0]-[xc(1,:) 0]),([xc(3,:) 0]-[xc(1,:) 0]))/...
        norm(xc(1,:)-xc(2,:))/norm(xc(1,:)-xc(3,:));
    
    if eta(3) < 0
        temp = rotIndx(elem,2);
        rotIndx(elem,2) = rotIndx(elem,3);
        rotIndx(elem,3) = temp;        
    end
end


end