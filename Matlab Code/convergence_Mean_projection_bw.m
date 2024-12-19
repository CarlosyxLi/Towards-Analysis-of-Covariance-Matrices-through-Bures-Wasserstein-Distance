function [out_mean,out_it,out_ssd] = convergence_Mean_projection_bw(X,eps)

% compute inital mean
    mean_int = mean(X,3);
    szx =size(X);
    num_matri = szx(3);
% Compute the sum sqaure of BW distance between mean_int and all 10 initial
% matrices
SSD = 0;
for m = 1:num_matri
    SSD = SSD + (compute_W_distance(mean_int, X(:,:,m)))^2;
end
       
% precision of algorithm
 k = 1;
 dist_m = 10;
    
 while dist_m > eps
     mean_k = mean_int;
     % projection to tangent space
     T_X = X;
     for i = 1:1:num_matri
          A = X(:,:,i);
          B = BW_log(mean_int,A);
          T_X(:,:,i) = B;
     end
        
     C = mean(T_X,3);
     % Project back to manifold
     mean_int = BW_exp(mean_k,C);
    
     %calculate the sum square of distance  
     ssd = 0;
     for m = 1:num_matri
         ssd = ssd + (compute_W_distance(mean_int, X(:,:,m)))^2;
     end
     SSD = [SSD, ssd];
    
     % Calculate BW distance
     dist_m = compute_W_distance(mean_k,mean_int);
        
     k = k + 1;
 end
out_it = i; % not necessary, computing the number of iteration
out_mean = mean_int; % mean for output
out_ssd = SSD;
end