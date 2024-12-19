function out_dis = compute_F_Distance(X,Y)

    out_dis = sqrt(trace((X-Y)*(X-Y)'));
end