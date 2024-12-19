function [out_mean,out_it] = compute_Mean_cheap2(X,eps,method)

sz = size(X);
center0 = mean(X,3);
maxdiff = ones(1,sz(3));
k = 0;
while any(maxdiff > eps)

    % BW
    if method == "BW"
        for i = 1:sz(3)
            % tangent point
            base = X(:,:,i);
            % projecting matrices to tangent space
            for j = 1:sz(3)
                Z(:,:,j) = BW_log(base,X(:,:,j));
            end
            % projecting arithmetic mean back to manifold
            C = mean(Z,3);
            X(:,:,i) = real(BW_exp(base,C));
        end
    
    % Riemannian
    else
        for i = 1:sz(3)
            % tangent point
            base = X(:,:,i);
            % projecting matrices to tangent space
            for j = 1:sz(3)
                Z(:,:,j) = Rie_log(base,X(:,:,j));
            end
            % projecting arithmetic mean back to manifold
            C = mean(Z,3);
            X(:,:,i) = real(Rie_exp(base,C));
        end

    end

    % update center
    center = mean(X,3);

    % abnormality?
    if mean(abs(center./center0),"all") > 100
        center = center0;
        center(:,:) = 100;
        k = -1;
        break
    end

    % convergence - maximal element difference between each matrix and center
    for i = 1:sz(3)
        maxdiff(1,i) = max(max(abs(X(:,:,i)-center)));
    end

    k = k+1;

end

out_mean = center;
out_it = k;

end