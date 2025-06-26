function [totalE, P_global, K_global] = ...
    magnetAssemble(shapeFunction, indx, UUR, magnetInfo, magnetConnection, ...
    totalE,P_global, K_global, ...
    targetElemID, XcmList, etaList, Ng1, Ng2,...
    dfctr, calP, calK)

global mu0

% numNodeElem =  size(indx,2);
% nDof =  size(xNode,1) * DIM;

xcmList = calCM(shapeFunction, indx, targetElemID,XcmList, etaList, UUR, dfctr);


for mag = 1:size(magnetConnection,1)
    mag1 = magnetConnection(mag,1);
    mag2 = magnetConnection(mag,2);
    xcm1 =  xcmList(mag1,:);
    xcm2 =  xcmList(mag2,:);
    m1 = magnetInfo{mag1,3};
    m2 = magnetInfo{mag2,3};
    d = norm(xcm1 - xcm2);

    E  = -mu0*m1*m2/4/pi/norm(xcm1 - xcm2)^3;
    totalE = totalE + E;

    if calP
        F1 = -3*E/d^2*(xcm1' - xcm2');
        F2 = -F1;

        P_global = P_global + Ng1'*F1 + Ng2'*F2;
    end
    
    if calK
    [K11, K12, K22] = d2EdVxcmdVxcm(E,xcm1,xcm2);

    K_global = K_global + Ng1'* K11*Ng1 + ...
        Ng2'* K22*Ng2 + Ng1'* K12*Ng2 + ...
        Ng2'* K12*Ng1;
    

    end
end

end

function [K11, K12, K22] = d2EdVxcmdVxcm(E,Vxcm1,Vxcm2)

K11 = eye(3);
K11 =  sparse(K11);
d = norm(Vxcm1 - Vxcm2);


K11 =  K11.*(15*E/d^4*(Vxcm1'-Vxcm2').^2 - 3*E/d^2);

K11(1,2) = 15*E/d^4*(Vxcm1(1) - Vxcm2(1))*(Vxcm1(2) - Vxcm2(2));
K11(2,1) =  K11(1,2);

K11(2,3) = 15*E/d^4*(Vxcm1(2) - Vxcm2(2))*(Vxcm1(3) - Vxcm2(3));
K11(3,2) = K11(2,3);

K11(3,1) = 15*E/d^4*(Vxcm1(3) - Vxcm2(3))*(Vxcm1(1) - Vxcm2(1));
K11(1,3) = K11(3,1);

K22 = K11;
K12 = -K11;

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
