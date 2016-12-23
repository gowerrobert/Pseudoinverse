function uniform_sym_rank(n,r)
    A = randn(n,n);   # randn(n,n)
    A = A'+A;
(U, S, V)= svd(A);
    A= U[:, 1: r]* diagm(S[1: r])* U[:, 1: r]';
    sol = randn(n,1);
    b = A*sol;
# now replace with least norm solution
    pinvA =  pinv(A); 
    sol = pinvA*b;
    title = string("uniform-random-sym-" , string(n) , "_r_" , string(r))
    return Prob(A,b,sol, pinvA, title)
end



#A = randn(n);   % randn(n,n)
#A = A+A';
#[U, S, ~]= svd(A);
#B= U(:, 1: r)* S(1: r, 1: r)* U(:, 1: r)';
#Prob.sol = randn(n,1);
#Prob.b = B*Prob.sol;
#% now replace with least norm solution
#Prob.sol = pinv(B)*Prob.b;
#Prob.A =B;
#Prob.title =[ 'uniform-random_n_'  num2str(n) '_r_' num2str(r)];