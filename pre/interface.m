function [indx, rotIndx] = interface(k,kt,pattern,N1, N2_list, xCoord, indx,rotIndx,ucID)

N2 = 0;
nElem = size(indx,1);
nRot = 0;

for n =  1: length(pattern)-1
    N2 = N2 + N2_list(n);
    
    %%% interface axial springs
    
    for p = 1:N1
        nElem  = nElem + 1;
        n1 = ucID(p,N2,4);
        n2 = ucID(p,N2+1,1);
        indx(nElem,1:2) =  [n1 n2];
        indx(nElem,3) = k;
        indx(nElem,4) = norm(xCoord(n1,1:2)-xCoord(n2,1:2));
    end
    
    %define rotational springs
    
    for p = 1: N1
%         nRot =  nRot + 1;
        n1 = ucID(p,N2,4);
        n2 = ucID(p,N2,3);
        n3 = ucID(p,N2+1,1);
        rotIndx= [rotIndx;n1 n2 n3 kt];
%         
%         nRot =  nRot + 1;
        n1 = ucID(p,N2+1,1);
        n2 = ucID(p,N2+1,2);
        n3 = ucID(p,N2,4);
        rotIndx= [rotIndx;n1 n2 n3 kt];
        
        if p ~= 1
%         nRot =  nRot + 1;
        n1 = ucID(p,N2,4);
        n2 = ucID(p-1,N2,3);
        n3 = ucID(p,N2+1,1);
        rotIndx= [rotIndx;n1 n2 n3 kt];
        
%         nRot =  nRot + 1;
        n1 = ucID(p,N2+1,1);
        n2 = ucID(p-1,N2+1,2);
        n3 = ucID(p,N2,4);
        rotIndx= [rotIndx;n1 n2 n3 kt];
        
        end       
                       
    end
end

end