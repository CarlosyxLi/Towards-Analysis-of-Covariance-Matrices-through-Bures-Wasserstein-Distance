function out_dis = compute_F_norm(X)

    out_dis = sqrt(trace(X*X'));
end