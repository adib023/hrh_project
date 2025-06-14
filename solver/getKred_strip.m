function Kred = getKred_strip(K,W,Wint,mu, right_rows, left_rows)

Imat = eye(size(K));

for p = right_rows
    Imat(p,p) = exp(i*mu);
end

for p = left_rows
    Imat(p,p) = exp(-i*mu);
end

Wmat =  Imat * W ;

Kred =  Wmat' * K * Wint;

end