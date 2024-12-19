function [out_mean,out_it,out_ssd] = convergence_Mean_inductive(X,eps,method)
%Convergence of Inductive mean for BW
% put matrices into array
szx =size(X);
num_matri = szx(3);
matrixtmp1 = X(:,:,1)*X(:,:,2);
matrixtmp1 = modify_negative_matrix(matrixtmp1);    
matrixtmp2 = X(:,:,2)*X(:,:,1);
matrixtmp2 = modify_negative_matrix(matrixtmp2);
    
% compute inital mean by geodesic of bw method
t=1/2;
mean_int =  ((1-t)^2)*X(:,:,1)+(t^2)*X(:,:,2)+t*(1-t)*((matrixtmp1)^(1/2)+(matrixtmp2)^(1/2));

% Compute the sum sqaure of BW distance between mean_int and all 10 initial
% matrices

SSD = 0;
for m = 1:num_matri
    if method == "BW"
        SSD = SSD + (compute_W_distance(mean_int, X(:,:,m)))^2;
    else
        SSD = SSD + (compute_Rie_distance(mean_int, X(:,:,m)))^2;
    end        
end

k = 2;
i = 1;
o = 3;
dist_m = 10;
    
% the iteration of means by geodesic of bw method
while dist_m > eps
    k = k + 1;
    mean_k = mean_int;
    mat = X(:,:,o);
        
    matrixtmp3 = mean_int*mat;
    matrixtmp3 = modify_negative_matrix(matrixtmp3);
    matrixtmp4 = mat*mean_int;
    matrixtmp4 = modify_negative_matrix(matrixtmp4);

    mean_int =  ((1-1/k)^2)*mean_int+((1/k)^2)*mat+(1/k)*(1-1/k)*((matrixtmp3)^(1/2)+(matrixtmp4)^(1/2));
    mean_k1 = mean_int;
    if method == "BW"
        dist_m = compute_W_distance(mean_k,mean_k1);
    else
        dist_m = compute_Rie_distance(mean_k,mean_k1);
    end
        
    %calculate the sum square of distance  
    ssd = 0;
    for m = 1:num_matri
        if method == "BW"
            ssd = ssd + (compute_W_distance(mean_k1, X(:,:,m)))^2;
        else
            ssd = ssd + (compute_Rie_distance(mean_k1, X(:,:,m)))^2;
        end
    end
    SSD = [SSD, ssd];
    % computing the number of iteration
    if o == num_matri 
        o = 1;
    else
        o = o+1;
    end
    i = i + 1;
end
out_it =i; % not necessary, computing the number of iteration
out_mean = mean_k1; % mean for output
out_ssd = SSD;
end

