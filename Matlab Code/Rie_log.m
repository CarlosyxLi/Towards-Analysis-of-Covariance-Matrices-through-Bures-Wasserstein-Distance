function out_dis = Rie_log(X,Y)
    out_dis = (X^(1/2))*logm((X^(-1/2))*Y*(X^(-1/2)))*(X^(1/2));
end