function gamma = compute_BW_geodesic_mean(X,Y)
%equation II.6 - BW geodesic mean of two matrices
    gamma = 1/4 * (X + Y + sqrtm(X*Y) + sqrtm(Y*X));
end