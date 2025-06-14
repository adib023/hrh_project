function rotIndx = fix_set_angle(xCoord,rotIndx)
for elem = 1:size(rotIndx,1)
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


for p = 1: size(rotIndx,1)
    n1 = rotIndx(p,1);
    n2 = rotIndx(p,2);
    n3 = rotIndx(p,3);
    
    xc = xCoord([n1 n2 n3],:);
    
    costh_0 =  dot((xc(2,:)-xc(1,:)),(xc(3,:)-xc(1,:)))/...
        norm(xc(1,:)-xc(2,:))/norm(xc(1,:)-xc(3,:));
    sinth_0 =  norm(cross(([xc(2,:) 0]-[xc(1,:) 0]),([xc(3,:) 0]-[xc(1,:) 0])))/...
        norm(xc(1,:)-xc(2,:))/norm(xc(1,:)-xc(3,:));
    th_0 =  atan2(sinth_0,costh_0);
    rotIndx(p,5) = th_0;
end

end