function out_dis = BW_log(X,Y)
   out_dis = real(real(sqrtm((X*Y)))+real(sqrtm((Y*X)))-2*X);
end