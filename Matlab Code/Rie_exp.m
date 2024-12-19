function out_dis = Rie_exp(A,X)
    out_dis = (A^(1/2))*expm((A^(-1/2))*X*(A^(-1/2)))*(A^(1/2));
end