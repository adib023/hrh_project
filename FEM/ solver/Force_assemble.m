function PF = Force_assemble(n_load, force_node, force_val, eq_num, PF)

% assemble PF

for i = 1: n_load
    
    node = force_node(i,1);
    dof = force_node(i,2);
    f = force_val(i);
    row = eq_num(node,dof);
    if row > 0
        PF(row) = PF(row) + f;
    end
end

return;