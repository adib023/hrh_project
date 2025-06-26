function Kel =  d2VxcmdUdU_part(shapeFunction, nGaussPts, gaussPts, gaussWeights ...
    ,Xel,Uel,Xcm,dEdVxcm, nNodeElem)

global DIM

Kel = zeros(nNodeElem*DIM,nNodeElem*DIM);

for p = 1:nGaussPts
    r =  gaussPts(p,1:3);
    w = gaussWeights(p);

    [N,DN] = shapeFunction(r);
    J = Xel*DN;
    detJ =  det(J);

    invJ = inv(J);
    B = (DN*invJ);
    NoI = kron(N',eye(3))';
    defF = Uel*B + eye(3);
    deterF =  det(defF);
    BoI = kron(B',eye(3));
    cofF = (deterF*inv(defF))';
    cofF = cofF(:); %%%%% 9 by 1
    delCofF =  cal_del_cofF_del_u( defF, BoI );

    Kel = Kel +  (NoI' * dEdVxcm) * (BoI' * cofF)' * detJ * w +...
        ((BoI' * cofF) * (dEdVxcm' * NoI))  * detJ * w + ...
        (dEdVxcm'*(Xcm + Uel*N')) * ( BoI' * delCofF ) * detJ * w;

end

end

function delCofF = cal_del_cofF_del_u( F, BoI )

FF = squareTensorOpt(inv(F)',inv(F));
cofF  = det(F) * (inv(F)');
FinvT = inv(F)';

T = transposer ;
delCofF = FinvT(:) * (cofF(:)'*BoI) - det(F)*(FF*T*BoI);
end

function [C2] = transposer()

C4 = zeros(3,3,3,3) ;
C4(1,1,1,1) = 1 ; C4(1,2,2,1) = 1 ; C4(1,3,3,1) = 1 ;
C4(2,1,1,2) = 1 ; C4(2,2,2,2) = 1 ; C4(2,3,3,2) = 1 ;
C4(3,1,1,3) = 1 ; C4(3,2,2,3) = 1 ; C4(3,3,3,3) = 1 ;

for p = 1:3
    for q = 1:3
        for r = 1:3
            for s = 1:3
                indx1 = 3*(q-1)+p ;
                indx2 = 3*(s-1)+r ;

                C2(indx1,indx2) = C4(p,q,r,s);
            end
        end
    end
end

end

function [C2] = squareTensorOpt(A,B)

for p = 1:3
    for q = 1:3
        for r = 1:3
            for s = 1:3
                %%%%%%% according to book, it should be
                C4(p,q,r,s) = A(p,r)*B(q,s);

                indx1 = 3*(q-1)+p ;
                indx2 = 3*(s-1)+r ;

                C2(indx1,indx2) = C4(p,q,r,s);
            end
        end
    end
end

end



