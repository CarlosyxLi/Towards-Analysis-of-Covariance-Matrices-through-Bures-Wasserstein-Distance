function out_dis = compute_Rie_distance(X,Y)
    lambda = real(eig(Y/X));
    tmp_dis = sum((log(lambda)).^2);
    out_dis = real(sqrt(tmp_dis));
end