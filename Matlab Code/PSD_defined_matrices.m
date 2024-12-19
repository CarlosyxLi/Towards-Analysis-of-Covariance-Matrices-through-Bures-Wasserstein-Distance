%num_dimen = 5
%numzero = 2

function out_sigma1 = PSD_defined_matrices(num_dimen,numzero,eps1)
    A=rand(num_dimen,(num_dimen-numzero));
    A=A*A';
    B=eps1*diag(rand(1,num_dimen));
    out_sigma1 = A+B;
end
    
