function out_dis = BW_exp(A,X)
[U,D] = eig(A);
Xu = U'*X*U;
W = D;
for i = 1:size(D,1)
    for j = 1:size(D,2)
        W(i,j) = 1/(D(i,i)+D(j,j));
    end
end
C = (W.*Xu)*D*(W.*Xu);
out_dis = A+X+U*(C)*U';
end