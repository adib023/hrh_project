function w = hexa_disp_anal(xCoord,indx,rotIndx,ka,kt,k)

% ka - axial stiffness
% kt - torsional stiffness
% k -  wavevector


K = zeros(4);

for elem = 1:3
    n1 = indx(elem,1);
    n2 = indx(elem,2);
    
    xv1 = xCoord(n1,:);
    xv2 = xCoord(n2,:);
    
    Ka = getAxialSpringK(xv1,xv2,ka);
    size(Ka)
    
    kDotxlp =  dot(k,xCoord(n2,:) - xCoord(1,:))
    
    Ka(1:2,3:4) = Ka(1:2,3:4)*exp(-i*kDotxlp);
    Ka(3:4,1:2) = Ka(3:4,1:2)*exp(i*kDotxlp);
    K = K + Ka;
end


for elem = 1:3
  nodes = rotIndx(elem,:);
  xc=xCoord(nodes,:);
  
  Kt =  getRotK(xc,kt);
  
  kDotxlp1 =  dot(k,xCoord(nodes(2),:) - xCoord(1,:));
  kDotxlp2 =  dot(k,xCoord(nodes(3),:) - xCoord(1,:)); 
  
  Kt(1:2,3:4) = Kt(1:2,3:4)*exp(-i*kDotxlp1);
  Kt(1:2,5:6) = Kt(1:2,5:6)*exp(-i*kDotxlp2);
  
  K(1:2,1:2) = K(1:2,1:2)+ Kt(1:2,1:2);
  K(1:2,3:4) = K(1:2,3:4)+ Kt(1:2,3:4) + Kt(1:2,5:6);
end



for elem = 1:3
  nodes = rotIndx(elem,:);
  xc=xCoord(nodes,:);
  
  Kt =  getRotK(xc,kt);
  
  kDotxlp1 =  dot(k,xCoord(nodes(2),:) - xCoord(1,:));
  kDotxlp2 =  dot(k,xCoord(nodes(3),:) - xCoord(1,:)); 
  
  Kt(1:2,3:4) = Kt(1:2,3:4)*exp(i*kDotxlp1);
  Kt(1:2,5:6) = Kt(1:2,5:6)*exp(i*kDotxlp2);
  
  K(3:4,3:4) = K(3:4,3:4)+ Kt(1:2,1:2);
  K(3:4,1:2) = K(3:4,1:2)+ Kt(1:2,3:4) + Kt(1:2,5:6);
end


w2 = eig(K);

w = sqrt(w2);
w = real(w);
w = sort(w);


end
