function uniform_mat_rank(m,n,r)
A = randn(m,n);   # randn(n,n)
(U, S, V)= svd(A);
A= U[:, 1: r]* diagm(S[1: r])* V[:, 1: r]';
sol = randn(n,1);
b = A*sol;
# now replace with least norm solution
pinvA =  pinv(A); 
sol = pinvA*b;
title = string("uniform-random" , string(m) , "X" , string(n) , "_r_" , string(r))
    return Prob(A,b,sol, pinvA, title)
end