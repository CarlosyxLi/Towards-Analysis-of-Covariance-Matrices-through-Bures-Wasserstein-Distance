function out_dis = compute_W_distance(X,Y)
        tmp_dis = trace(X + Y)-2*trace(real(sqrtm((X*Y)))); 
        out_dis = real(sqrt(tmp_dis));
end