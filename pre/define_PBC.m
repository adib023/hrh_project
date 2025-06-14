function [axialpbcIndx, rotPBC] = define_PBC(N1,N2,xCoord,indx,rotIndx,ucID)
%%%%% define periodic boundary element

nAxPBC = 0;

for p = 1: N2
    nAxPBC = nAxPBC + 1;
    n1 = ucID(1,p,1);
    n2 = ucID(N1,p,2);
    n3 = ucID(N1,p,1);
    n4 = ucID(1,p,2);
    [k,l_0] = checkL(indx,n1,n4);
    axialpbcIndx(nAxPBC,1) = n1;
    axialpbcIndx(nAxPBC,2) = n2;
    axialpbcIndx(nAxPBC,3) = k;
    axialpbcIndx(nAxPBC,4:5) = (xCoord(n3,:)-xCoord(n1,:))/(N1-1) * N1;
    axialpbcIndx(nAxPBC,6:7) = [0 0];
    axialpbcIndx(nAxPBC,8) = l_0;
    
    nAxPBC = nAxPBC + 1;
    n1 = ucID(1,p,4);
    n2 = ucID(N1,p,3);
    n3 = ucID(N1,p,4);
    n4 = ucID(1,p,3);
    [k,l_0] = checkL(indx,n4,n1);
    axialpbcIndx(nAxPBC,1) = n1;
    axialpbcIndx(nAxPBC,2) = n2;
    axialpbcIndx(nAxPBC,3) = k;
    axialpbcIndx(nAxPBC,4:5) = (xCoord(n3,:)-xCoord(n1,:))/(N1-1) * N1;
    axialpbcIndx(nAxPBC,6:7)= [0 0];
    axialpbcIndx(nAxPBC,8) = l_0;
    
end






numRotPBC =  0;
for p = 1:N2
    numRotPBC = numRotPBC + 1;
    n1 = ucID(N1,p,2);
    n2 = ucID(N1,p,1);
    n3 = ucID(1,p,1);
    n_check = [ucID(1,p,2),ucID(1,p,1),ucID(2,p,1)];
    l1 = [0 0];
    l2 = [0 0];
    l3 = (xCoord(n2,:)-xCoord(n3,:))/(N1-1) * N1;
    [kt,th_0] = checkTheta(rotIndx, n_check);
    rotPBC(numRotPBC,:) = [n1 n2 n3 kt l1 l2 l3 th_0];
    
    numRotPBC = numRotPBC + 1;
    n1 = ucID(N1,p,2);
    n2 = ucID(N1,p,3);
    n3 = ucID(1,p,1);
    n_check = [ucID(1,p,2),ucID(1,p,3),ucID(2,p,1)];
    l1 = [0 0];
    l2 = [0 0];
    l3 = (xCoord(n1,:)-xCoord(ucID(1,p,2),:))/(N1-1) * N1;
    [kt,th_0] = checkTheta(rotIndx, n_check);
    rotPBC(numRotPBC,:) = [n1 n2 n3 kt l1 l2 l3 th_0];
    
    
    numRotPBC = numRotPBC + 1;
    n1 = ucID(N1,p,3);
    n2 = ucID(N1,p,2);
    n3 = ucID(1,p,4);
    n_check = [ucID(1,p,3),ucID(1,p,2),ucID(2,p,4)];
    l1 = [0 0];
    l2 = [0 0];
    l3 = (xCoord(n2,:)-xCoord(ucID(1,p,2),:))/(N1-1) * N1;
    [kt,th_0] = checkTheta(rotIndx, n_check);
    rotPBC(numRotPBC,:) = [n1 n2 n3 kt l1 l2 l3 th_0];
    
    numRotPBC = numRotPBC + 1;
    n1 = ucID(N1,p,3);
    n2 = ucID(N1,p,4);
    n3 = ucID(1,p,4);
    n_check = [ucID(1,p,3),ucID(1,p,4),ucID(2,p,4)];
    l1 = [0 0];
    l2 = [0 0];
    l3 = (xCoord(n1,:)-xCoord(ucID(1,p,3),:))/(N1-1) * N1;
    [kt,th_0] = checkTheta(rotIndx, n_check);
    rotPBC(numRotPBC,:) = [n1 n2 n3 kt l1 l2 l3 th_0];
    
    
    numRotPBC = numRotPBC + 1;
    n1 = ucID(1,p,1);
    n2 = ucID(1,p,2);
    n3 = ucID(N1,p,2);
    n_check = [ucID(2,p,1),ucID(2,p,2),ucID(1,p,2)];
    l1 = (xCoord(n3,:)-xCoord(n2,:))/(N1-1) * N1;
    l2 = (xCoord(n3,:)-xCoord(n2,:))/(N1-1) * N1;
    l3 = [0 0];
    [kt,th_0] = checkTheta(rotIndx, n_check);
    rotPBC(numRotPBC,:) = [n1 n2 n3 kt l1 l2 l3 th_0];
    
    if p ~= 1
    numRotPBC = numRotPBC + 1;
    n1 = ucID(1,p,1);
    n2 = ucID(1,p-1,4);
    n3 = ucID(N1,p,2);
    n_check = [ucID(2,p,1),ucID(1,p,2),ucID(2,p-1,4)];
    l1 = (xCoord(n3,:)-xCoord(ucID(1,p,2),:))/(N1-1) * N1;
    l2 = (xCoord(n3,:)-xCoord(ucID(1,p,2),:))/(N1-1) * N1;
    l3 = [0 0];
    [kt,th_0] = checkTheta(rotIndx, n_check);
    rotPBC(numRotPBC,:) = [n1 n2 n3 kt l1 l2 l3 th_0];
    end
    
    
    numRotPBC = numRotPBC + 1;
    n1 = ucID(1,p,4);
    n2 = ucID(1,p,3);
    n3 = ucID(N1,p,3);
    n_check = [ucID(2,p,4),ucID(2,p,3),ucID(1,p,3)];
    l1 = (xCoord(n3,:)-xCoord(n2,:))/(N1-1) * N1;
    l2 = (xCoord(n3,:)-xCoord(n2,:))/(N1-1) * N1;
    l3 = [0 0];
    [kt,th_0] = checkTheta(rotIndx, n_check);
    rotPBC(numRotPBC,:) = [n1 n2 n3 kt l1 l2 l3 th_0];
    
    if p ~= N2
    numRotPBC = numRotPBC + 1;
    n1 = ucID(1,p,4);
    n2 = ucID(1,p+1,1);
    n3 = ucID(N1,p,3);
    n_check = [ucID(2,p,4),ucID(1,p,3),ucID(2,p+1,1)];
    l1 = (xCoord(n3,:)-xCoord(ucID(1,p,3),:))/(N1-1) * N1;
    l2 = (xCoord(n3,:)-xCoord(ucID(1,p,3),:))/(N1-1) * N1;
    l3 = [0 0];
    [kt,th_0] = checkTheta(rotIndx, n_check);
    rotPBC(numRotPBC,:) = [n1 n2 n3 kt l1 l2 l3 th_0];
    end
end




end

function [k,L] = checkL(indx,n1,n2)
for p = 1 : size(indx,1)
    if indx(p,1) == n1 & indx(p,2) == n2
        k = indx(p,3);
        L = indx(p,4);
    end
end
end


function [kt,theta] = checkTheta(rotIndx, nodes)
for p = 1:size(rotIndx,1)
    if isequal(sort(rotIndx(p,1:3)),sort(nodes))
        kt = rotIndx(p,4);
        theta =  rotIndx(p,5);
    end
end
end

