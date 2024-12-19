% Number of matrices
matcount = (10:10:300); 
% Dimension setting
dimcount = (5:15:100);
%stop criteria
eps = 0.001;


%Generate save matrix
BW_induc = zeros(length(dimcount),length(matcount));
BW_proje = zeros(length(dimcount),length(matcount));
BW_cheap = zeros(length(dimcount),length(matcount));
Rie_induc = zeros(length(dimcount),length(matcount));
Rie_proje = zeros(length(dimcount),length(matcount));
Rie_cheap = zeros(length(dimcount),length(matcount));

%Start calcualtion
%%
for n = 6:length(dimcount)
    for numArrays = 1:length(matcount)
        
        %generate matrices
        X = random_positivedefined_matrices(dimcount(n));
        for i = 2:matcount(numArrays)
            X(:,:,i) = random_positivedefined_matrices(dimcount(n));
        end
        
        %compute 3 Rie means
        f = @() compute_Mean_inductive(X,eps,"Rie");
        Rie_induc(n,numArrays) = timeit(f);

        f = @() compute_Mean_projection(X,eps,"Rie");
        Rie_proje(n,numArrays) = timeit(f);  

        f = @() compute_Mean_cheap(X,eps,"Rie");
        Rie_cheap(n,numArrays) = timeit(f);

        %compute 3 BW means
        f = @() compute_Mean_inductive(X,eps,"BW");
        BW_induc(n,numArrays) = timeit(f);

        f = @() compute_Mean_projection(X,eps,"BW");
        BW_proje(n,numArrays) = timeit(f);  

        f = @() compute_Mean_cheap(X,eps,"BW");
        BW_cheap(n,numArrays) = timeit(f);
    end
end




