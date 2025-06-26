function [TotalEnergy,P_global,K_global,detVol,elim_eqn] = ...
    elasticAssemble (ndim,n_elem, coords, elem_node, mat_stiff,elemMatIndx,...
                periodic_node, fixNode, eq_num, nDof, UUR)
TotalEnergy = 0 ; 
detVol = 1 ; 
for elem_num =  1:n_elem %5448
    matID =  elemMatIndx(elem_num);
    elem_stiff = mat_stiff(matID,:);

    [Kel, Pel,ElemEnergy,detVol]   =  element('STIFFNESSMATRIX',ndim, elem_num, coords, elem_node, ...
                    elem_stiff, UUR);

    [K_global,P_global] = element_stiffness_assemble(ndim,elem_node,eq_num,...
                        elem_num,Kel,Pel,K_global,P_global);
                    
	if (detVol == 0)
		break; 
	end

	TotalEnergy = TotalEnergy + ElemEnergy ;
end

%%%%%%%% Periodic Boundary condition reduction: 
%%%%%%doing row sum
elim_eqn = [];
for  p = 1:size(periodic_node,1)
    master = periodic_node(p,1);
    follower = periodic_node(p,2);
    for d = 1:ndim
        follow_eqn = eq_num(follower,d);
        master_eqn = eq_num(master,d);
        K_global(master_eqn,:) =  K_global(master_eqn,:)+ K_global(follow_eqn,:);
        P_global(master_eqn,:) =  P_global(master_eqn,:)+ P_global(follow_eqn,:);
        elim_eqn = [elim_eqn follow_eqn];
    end
end

%%%% doing column sum
for  p = 1:size(periodic_node,1)
    master = periodic_node(p,1);
    follower = periodic_node(p,2);
    for d = 1:ndim
        follow_eqn = eq_num(follower,d);
        master_eqn = eq_num(master,d);
        K_global(:,master_eqn) =  K_global(:,master_eqn)+ K_global(:,follow_eqn);
    end
end


%%%%% apply fixed condition in K_global and P_global

elim_eqn = [elim_eqn eq_num(fixNode,1:ndim)];

%%%% eliminate all equation for z direction

% for n= 1:size(coords,1)
%     
%     eq = eq_num(n,3);
% 
%     if ~ismember(eq,elim_eqn)
%         
%         elim_eqn = [elim_eqn eq];
%     end
% end
    
K_global(elim_eqn,:)= [];
K_global(:,elim_eqn)= [];
P_global(elim_eqn) = [];

end
