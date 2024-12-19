function [out_mean,out_dist,out_it] = compute_Mean_projection2(X,eps1,eps2,seq,method)

sz = size(X);
% initial mean
S(:,:,1) = mean(X,3);

dist = eps+1; % initializing distance array
fluc = eps+1; % initializing fluctuation array
k = 2;
while dist(length(dist)) > eps1 & sum(fluc > eps2)
    % stopping criteria 1: distance between two successive means is smaller than eps1
    % stopping criteria 2: fluctuations of the last seq number of distances are all smaller than eps2

    % BW
    if method == "BW"
        % projecting matrices to tangent space
        for i = 1:sz(3)
            % updating matrices
            Z(:,:,i) = BW_log(S(:,:,k-1),X(:,:,i));
        end
        % projecting new mean back to manifold
        S(:,:,k) = BW_exp(S(:,:,k-1),mean(Z,3));
        % distance between initial mean and latest mean 
        dist(k-1) = compute_W_distance(S(:,:,k-1),S(:,:,k));

    % Riemannian
    else
        % projecting matrices to tangent space
        for i = 1:sz(3)
            % updating matrices
            Z(:,:,i) = Rie_log(S(:,:,k-1),X(:,:,i));
        end
        % projecting new mean back to manifold
        S(:,:,k) = Rie_exp(S(:,:,k-1),mean(Z,3));
        % distance between initial mean and latest mean 
        dist(k-1) = compute_Rie_distance(S(:,:,k-1),S(:,:,k));

    end

    % fluctuation
    if length(dist) > seq % only computing when there are more than seq number of terms
        % compute fluctuation of the last seq number of distances
        fluc = dist(end-seq:end-1) - dist(end-seq+1:end);
        fluc = abs(fluc);
    end

    k = k+1;

end

out_mean = S(:,:,end);
out_dist = dist;
out_it = k-2;

end