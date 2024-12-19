function gamma = compute_Rie_geodesic_mean(X,Y)
%equation I.4 - Riemannian geodesic mean of two matrices
    gamma = sqrtm(X) * sqrtm(sqrtm(inv(X)) * Y * sqrtm(inv(X))) * sqrtm(X);
end