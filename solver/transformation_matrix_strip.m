%%% determine transformation matrix
function [W,Wint, left_rows, right_rows] = ... 
    transformation_matrix_strip(nNode, ...
    strip_fix,stripNode, rightNode, leftNode)


W =  zeros(nNode*2, (length(stripNode)-length(strip_fix))*2);
Wint = W ; 

col =  0;
left_rows = [];
right_rows = [];

for num = 1:length(stripNode)
    n =  stripNode(num)
    if (ismember(n,strip_fix)== false)
        id = find(stripNode == n);
        l =  leftNode(id)
        r =  rightNode(id)
    for m = 1:2 
            col = col + 1;
            row_s = 2*(n-1) + m;
            row_l = 2*(l-1) + m;
            row_r = 2*(r-1) + m;

            W(row_s,col) =  1;
            W(row_l,col) =  1;
            W(row_r,col) =  1;

            Wint(row_s , col) = 1 ; 

            left_rows = [left_rows row_l];
            right_rows = [right_rows row_r];
    end
    end
end
end