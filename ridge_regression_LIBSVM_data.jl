function ridge_regression_LIBSVM_data(dataset)
    X,y = loadDataset(dataset);
    A = X'*X;
#    sA = size(A);
#    sol = randn(sA[1],1);
#    b = A*sol;
#    pinvA =  pinv(A); 
#    sol = pinvA*b;
    title = "$(dataset)-ridge";
    return Prob(A,[],[], [], title)
end