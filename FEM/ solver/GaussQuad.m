function [nGaussPts,gaussPts,gaussWts] = GaussQuad(dim,elemType)

% sType : simplex type 
%	1: tetrahedron
%	2: hexahedron


if (dim  == 3) && (strcmp(elemType,'TET1') == 1)

	nGaussPts	= 4 ;

	gaussPts(1,1) = 0.5854101966249685;
	gaussPts(1,2) = 0.1381966011250105;
	gaussPts(1,3) = 0.1381966011250105;

	gaussPts(2,1) = 0.1381966011250105;
	gaussPts(2,2) = 0.5854101966249685;
	gaussPts(2,3) = 0.1381966011250105;

	gaussPts(3,1) = 0.1381966011250105;
	gaussPts(3,2) = 0.1381966011250105;
	gaussPts(3,3) = 0.5854101966249685;

	gaussPts(4,1) = 0.1381966011250105;
	gaussPts(4,2) = 0.1381966011250105;
	gaussPts(4,3) = 0.1381966011250105;

	gaussWts(1:4) = 0.25/6.0 ;

elseif (dim == 3) && (strcmp(elemType,'TET2') == 1)

	nGaussPts	=	11 ;

	gaussPts(1,1) = 0.2500;
	gaussPts(1,2) = 0.2500;
	gaussPts(1,3) = 0.2500;

	gaussPts(2,1) = 0.714285714285714285e-1;
	gaussPts(2,2) = 0.714285714285714285e-1;
	gaussPts(2,3) = 0.714285714285714285e-1;

	gaussPts(3,1) = 0.785714285714285714;
	gaussPts(3,2) = 0.714285714285714285e-1;
	gaussPts(3,3) = 0.714285714285714285e-1;

	gaussPts(4,1) = 0.714285714285714285e-1;
	gaussPts(4,2) = 0.785714285714285714;
	gaussPts(4,3) = 0.714285714285714285e-1;

	gaussPts(5,1) = 0.714285714285714285e-1;
	gaussPts(5,2) = 0.714285714285714285e-1;
	gaussPts(5,3) = 0.785714285714285714;

	gaussPts(6,1) = 0.399403576166799219;
	gaussPts(6,2) = 0.399403576166799219;
	gaussPts(6,3) = 0.100596423833200785;

	gaussPts(7,1) = 0.399403576166799219;
	gaussPts(7,2) = 0.100596423833200785;
	gaussPts(7,3) = 0.399403576166799219;

	gaussPts(8,1) = 0.399403576166799219;
	gaussPts(8,2) = 0.100596423833200785;
	gaussPts(8,3) = 0.100596423833200785;

	gaussPts(9,1) = 0.100596423833200785;
	gaussPts(9,2) = 0.399403576166799219;
	gaussPts(9,3) = 0.399403576166799219;

	gaussPts(10,1) = 0.100596423833200785;
	gaussPts(10,2) = 0.399403576166799219;
	gaussPts(10,3) = 0.100596423833200785;

	gaussPts(11,1) = 0.100596423833200785;
	gaussPts(11,2) = 0.100596423833200785;
	gaussPts(11,3) = 0.399403576166799219;

	gaussWts(1)  = -0.131555555555555550e-1;
	gaussWts(2)  =  0.762222222222222222e-2;
	gaussWts(3)  =  0.762222222222222222e-2;
	gaussWts(4)  =  0.762222222222222222e-2;
	gaussWts(5)  =  0.762222222222222222e-2;
	gaussWts(6)  =  0.248888888888888880e-1;
	gaussWts(7)  =  0.248888888888888880e-1;
	gaussWts(8)  =  0.248888888888888880e-1;
	gaussWts(9)  =  0.248888888888888880e-1;
	gaussWts(10) =  0.248888888888888880e-1;
	gaussWts(11) =  0.248888888888888880e-1;

elseif (dim == 3) && (strcmp(elemType,'HEX1') == 1)

	nGaussPts	= 8 ;

	gaussPts	= [-1/sqrt(3) , -1/sqrt(3) , -1/sqrt(3) ;  
			-1/sqrt(3) , -1/sqrt(3) , 1/sqrt(3) ;  
			-1/sqrt(3) ,  1/sqrt(3) ,  1/sqrt(3) ;  
			-1/sqrt(3) ,  1/sqrt(3) , -1/sqrt(3) ;  
			 1/sqrt(3) , -1/sqrt(3) , -1/sqrt(3) ;  
			 1/sqrt(3) , -1/sqrt(3) , 1/sqrt(3) ;  
			 1/sqrt(3) ,  1/sqrt(3) ,  1/sqrt(3) ;  
			 1/sqrt(3) ,  1/sqrt(3) , -1/sqrt(3)];  

	gaussWts(1:8) = 1.0/8.0 ;

elseif (dim == 3) && (strcmp(elemType,'HEX2') == 1)

	nGaussPts	=	27 ; 

	gp(1) = 0 ; gp(2) = -sqrt(3.0/5.0) ; gp(3) = -gp(2) ;  
	gwt(1) = 4/9.0 ; gwt(2) = 5/18.0 ; gwt(3) = 5/18.0 ; 

	gaussPts	= [	gp(1) , gp(1) , gp(1) ; 
				gp(1) , gp(1) , gp(2) ; 
				gp(1) , gp(1) , gp(3) ; 
				gp(1) , gp(2) , gp(1) ; 
				gp(1) , gp(2) , gp(2) ; 
				gp(1) , gp(2) , gp(3) ; 
				gp(1) , gp(3) , gp(1) ; 
				gp(1) , gp(3) , gp(2) ; 
				gp(1) , gp(3) , gp(3) ; 
				gp(2) , gp(1) , gp(1) ; 
				gp(2) , gp(1) , gp(2) ; 
				gp(2) , gp(1) , gp(3) ; 
				gp(2) , gp(2) , gp(1) ; 
				gp(2) , gp(2) , gp(2) ; 
				gp(2) , gp(2) , gp(3) ; 
				gp(2) , gp(3) , gp(1) ; 
				gp(2) , gp(3) , gp(2) ; 
				gp(2) , gp(3) , gp(3) ; 
				gp(3) , gp(1) , gp(1) ; 
				gp(3) , gp(1) , gp(2) ; 
				gp(3) , gp(1) , gp(3) ; 
				gp(3) , gp(2) , gp(1) ; 
				gp(3) , gp(2) , gp(2) ; 
				gp(3) , gp(2) , gp(3) ; 
				gp(3) , gp(3) , gp(1) ; 
				gp(3) , gp(3) , gp(2) ; 
				gp(3) , gp(3) , gp(3) ] ; 

	gaussWts = 8*[ gwt(1)*gwt(1)*gwt(1) , 
			 gwt(1)*gwt(1)*gwt(2) , 
			 gwt(1)*gwt(1)*gwt(3) , 
			 gwt(1)*gwt(2)*gwt(1) , 
			 gwt(1)*gwt(2)*gwt(2) , 
			 gwt(1)*gwt(2)*gwt(3) , 
			 gwt(1)*gwt(3)*gwt(1) , 
			 gwt(1)*gwt(3)*gwt(2) , 
			 gwt(1)*gwt(3)*gwt(3) , 
			 gwt(2)*gwt(1)*gwt(1) , 
			 gwt(2)*gwt(1)*gwt(2) , 
			 gwt(2)*gwt(1)*gwt(3) , 
			 gwt(2)*gwt(2)*gwt(1) , 
			 gwt(2)*gwt(2)*gwt(2) , 
			 gwt(2)*gwt(2)*gwt(3) , 
			 gwt(2)*gwt(3)*gwt(1) , 
			 gwt(2)*gwt(3)*gwt(2) , 
			 gwt(2)*gwt(3)*gwt(3) , 
			 gwt(3)*gwt(1)*gwt(1) , 
			 gwt(3)*gwt(1)*gwt(2) , 
			 gwt(3)*gwt(1)*gwt(3) , 
			 gwt(3)*gwt(2)*gwt(1) , 
			 gwt(3)*gwt(2)*gwt(2) , 
			 gwt(3)*gwt(2)*gwt(3) , 
			 gwt(3)*gwt(3)*gwt(1) , 
			 gwt(3)*gwt(3)*gwt(2) , 
			 gwt(3)*gwt(3)*gwt(3)];

else 
	error('quadrature not yet implemented') ; 
end

end
