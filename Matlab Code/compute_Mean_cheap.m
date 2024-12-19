function [out_mean,out_it] = compute_Mean_cheap(X,eps,method)

    szx =size(X);
    num_matri = szx(3);    
    meanarray_k = zeros(szx(1),szx(2),num_matri);
    mean_int = X;
    k = 0;
    dist_m = 10;

    while dist_m > eps % the maximum elements' difference among means

        if method == "BW"
        % compute means by every matrix's tangent space
            for i = 1:num_matri
                
                base = mean_int(:,:,i);
                B = mean_int;
                % projection to tangent space
                for j = 1:num_matri
                    A = mean_int(:,:,j);
                    B(:,:,j) = BW_log(base,A);
                end
                C = mean(B,3);
                meanarray_k(:,:,i) = real(BW_exp(base,C));
            end
        else
        % compute means by every matrix's tangent space
            for i = 1:num_matri
                
                base = mean_int(:,:,i);
                B = mean_int;
                % projection to tangent space
                for j = 1:num_matri
                    A = mean_int(:,:,j);
                    B(:,:,j) = Rie_log(base,A);
                end
                C = mean(B,3);
                meanarray_k(:,:,i) = Rie_exp(base,C);
            end
        end
        
         % determine convergence by the maximum element's difference
        center =  mean(real(meanarray_k),3);
        
        if method == "BW"
            dist_m = compute_W_distance(center,mean(real(mean_int),3));
        else
            dist_m = compute_Rie_distance(center,mean(real(mean_int),3));
        end
        
        %for i = 1:num_matri
         %   matabs = abs(meanarray_k(:,:,i)-center);
         %   matmax = max(matabs(:));
         %   dista(i,1)= matmax;
        %end

        mean_int =  real(meanarray_k);
        k = k +1;

    end
   out_it = k; % not necessary, computing the number of iteration
   out_mean = center; % mean for output

end