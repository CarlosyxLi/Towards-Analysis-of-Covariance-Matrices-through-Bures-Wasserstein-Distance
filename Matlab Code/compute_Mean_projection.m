function [out_mean,out_it] = compute_Mean_projection(X,eps,method)


% compute inital mean
    mean_int = mean(X,3);
    szx =size(X);
    num_matri = szx(3);
    
% precision of algorithm
    
    k = 1;
    dist_m = 10;
    
    while dist_m > eps
       
       mean_k = mean_int; %mean_k = S0

       if method == "BW"
           % projection to tangent space
           T_X = X; %input
           for i = 1:1:num_matri
                A = X(:,:,i);
                B = BW_log(mean_int,A); %project Ai onto tangent space at S0 and get Bi
                T_X(:,:,i) = B;
            end
            
            C = mean(T_X,3); %mean of Bis, C
            % Project back to manifold
            mean_int = BW_exp(mean_k,C); %project C onto tangent space at S0 and get mean_int = S1
            
            % Calculate distance
            dist_m = compute_W_distance(mean_k,mean_int); %distance between S0 and S1
       else
           % projection to tangent space
           T_X = X;
           for i = 1:1:num_matri
                A = X(:,:,i);
                B = Rie_log(mean_int,A);
                T_X(:,:,i) = B;
            end
            
            C = mean(T_X,3);
            % Project back to manifold
            mean_int = Rie_exp(mean_k,C);
            
            % Calculate distance
            dist_m = compute_Rie_distance(mean_k,mean_int);           
       end
       k = k + 1;
    
   end
   out_it = k; % not necessary, computing the number of iteration
   out_mean = mean_int; % mean for output
   
end