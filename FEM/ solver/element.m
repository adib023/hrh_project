function  varargout   =   element(shapeFunction, ...
    nGaussPts, gaussPts, gaussWeights,...
    ndim, elem_num, coords, elem_node, ...
    elem_stiff,rho, UUR, calP, calK, calM)


% [shapeFunction, nGaussPts, gaussWeights] =  getShapeFunction_and_gaussInfo(elemType);
nNodeElem = size(elem_node,2);


KEL = zeros(ndim*nNodeElem, ndim*nNodeElem);
PEL = zeros(ndim*nNodeElem, 1);
MEL =  zeros(ndim*nNodeElem, ndim*nNodeElem);


Xel = [coords(elem_node(elem_num,1:nNodeElem),1),  coords(elem_node(elem_num,1:nNodeElem),2), ...
    coords(elem_node(elem_num,1:nNodeElem),3)]';
Uel = [UUR(elem_node(elem_num,1:nNodeElem),1),  UUR(elem_node(elem_num,1:nNodeElem),2), ...
    UUR(elem_node(elem_num,1:nNodeElem),3) ]';


ElemEner = 0 ;
detVol = 1 ;

for p =1: nGaussPts

    r = gaussPts(p,1:ndim);
    w = gaussWeights(p) ;

    [N, DN] = shapeFunction(r);
    J = Xel*DN;

    detJ = det(J);

    if detJ < 0
        %Xel'
        %J
        %p
        elem_node(elem_num,:)
        elemID = elem_num
        %'determinant of J is negative!!!!'
    end

    invJ = inv(J);
    B = (DN*invJ);
    NoI = kron(N,eye(3));


    BoI = zeros(ndim*ndim,ndim*nNodeElem) ;

    BoI = kron(B',eye(3));


    defF = Uel*B + eye(ndim) ;
    deterF = det(defF) ;

    if (deterF < 0)
        % 'Negative detF'
        detVol = 0  ;
    end

    if ~calM
    [Ener , Tvec , Cmat] = Constitutive (elem_stiff,defF );
    ElemEner = ElemEner + detJ*w*Ener;
    end
    
    if calK
        KEL  = KEL + (BoI'*Cmat*BoI)*detJ*w;
	
    end
    if calP
        PEL  = PEL + (BoI'*Tvec)*detJ*w;
    end

    if calM
        MEL =  MEL + rho *( NoI'*NoI ) * detJ * w;
    end

end


if calM
    varargout{1} =  MEL;
else
    varargout{1} = KEL;
    varargout{2} = PEL;
    varargout{3} = ElemEner;
    varargout{4} = detVol;
end
end
