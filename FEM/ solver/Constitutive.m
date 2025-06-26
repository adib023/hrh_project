function [ener,Tstress,Cmat1] = Constitutive(elem_stiff,F)

[ener,T,Cmat1] = compNeoHookean(elem_stiff,F);

Tstress = T(:);


end

function [ener,T,Cmat1] = compNeoHookean(elem_stiff,F)

lamb = elem_stiff(1) ; mu = elem_stiff(2) ;

C = F'*F ;

ener = lamb/8.0*log(det(C))*log(det(C)) + ...
      0.5*mu*(C(1,1)+C(2,2)+C(3,3)-3-log(det(C)));

S = 0.5*lamb*log(det(C))*inv(C) + mu*(eye(3) - inv(C)) ;  
T = F*S ; 

C4mat2 = lamb * circ4(inv(C),inv(C)) + (2*mu-lamb*log(det(C)))*symsq4(inv(C),inv(C)) ;

delt = eye(3) ; 

for p = 1:3
	for q = 1:3
		for r = 1:3
			for s = 1:3
		
C4mat1(p,q,r,s) = delt(p,r)*S(q,s) ;  

		for t = 1:3
			for u = 1:3
C4mat1(p,q,r,s) = C4mat1(p,q,r,s) + F(p,t)*C4mat2(t,q,u,s)*F(r,u) ;
			end
		end
			end
		end
	end
end

Cmat1 = convert2Matrix(C4mat1) ;

end

function [C] = circ4(A,B)
  
  for p = 1:3
    for q = 1:3
      for r = 1:3
        for s = 1:3
          C(p,q,r,s) = A(p,q)*B(r,s);
        end
      end
    end
  end

end

function [C2] = convert2Matrix(C4)

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

function linearElasticIsotropic(elem_stiff) 

% linear elastic and isotropic: 

lamb = elem_stiff(1) ; 
mu   = elem_stiff(2) ;  

Dmat = zeros(6) ;

Dmat(1:3,1:3) = lamb*ones(3) + 2*mu*eye(3) ;  

Dmat(4:6,4:6) = mu*eye(3);

end

function [C] = symsq4(A,B)

for p = 1:3
	for q = 1:3
		for r = 1:3
			for s = 1:3 

C(p,q,r,s) = 0.25*(A(p,r)*B(q,s)+A(p,s)*B(q,r)+A(q,r)*B(p,s)+A(q,s)*B(p,r));
			
			end
		end
	end
end

end




