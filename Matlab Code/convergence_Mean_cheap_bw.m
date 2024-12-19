function [out_mean,out_it,out_ssd] = convergence_Mean_cheap_bw(X,eps)
    szx =size(X);
    num_matri = szx(3);    
    meanarray_k = zeros(szx(1),szx(2),num_matri);
    mean_int = X;
    k = 0;
    dista = ones(num_matri,1);
    SSD = 0;

    while all(dista(:) > eps)% the maximum elements' difference among means

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
            meanarray_k(:,:,i) = BW_exp(base,C);
        end
        
         % determine convergence by the maximum element's difference
        center =  mean(meanarray_k,3);
        ssd = 0;
        for i = 1:num_matri
            matabs = abs(meanarray_k(:,:,i)-center);
            matmax = max(matabs(:));
            dista(i,1)= matmax;
            %Calculate the sum square distance
            ssd = ssd + (compute_W_distance(center, X(:,:,i)))^2;
        end

        mean_int =  meanarray_k;
        SSD = [SSD, ssd];
        k = k +1;

    end
    
   out_it = k; % not necessary, computing the number of iteration
   out_mean = mean_int; % mean for output
   out_ssd = SSD(2:end);



end