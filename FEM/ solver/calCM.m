function xcmList = calCM(shapeFunction,indx, targetElemID, XcmList, etaList, UUR, dfctr)

global a1 a2 lr_pbc bt_pbc DIM

Xcm1 = XcmList(1,:);
Xcm2 = XcmList(2,:);

elem =  targetElemID(1);
elemNode =  indx(elem,:);
Uel =  UUR(elemNode,:)';

[N_cm, ~] =  shapeFunction(etaList(1,:));
NoI_cm =  kron(N_cm,eye(DIM));

xcm1 =  Xcm1 + (NoI_cm * Uel(:))';

elem =  targetElemID(2);
elemNode =  indx(elem,:);
Uel =  UUR(elemNode,:)';

[N_cm, ~] =  shapeFunction(etaList(2,:));
NoI_cm =  kron(N_cm,eye(DIM));


xcm2 =  Xcm2 + (NoI_cm * Uel(:))';

xcm3 = xcm1 + a1 + dfctr * bt_pbc;
xcm4 = xcm2 + a2 + dfctr * lr_pbc;

xcmList =  [xcm1; xcm2; xcm3; xcm4];


end
