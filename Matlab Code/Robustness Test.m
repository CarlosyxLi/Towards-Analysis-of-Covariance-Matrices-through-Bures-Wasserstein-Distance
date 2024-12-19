clear all;

%% parameters
e1 = 10^(-3); % 1st PSD calibration
e2 = 10^(-3); % strength of perturbation
e3 = 10^(-4); % 2nd PSD calibration (due to computational error)

eps = 0.001; % stopping criteria for inductive mean
eps1 = 0.001; % stopping criteria for BW projection mean
eps2 = 0.00001; % stopping criteria for Rie projection mean
seq = 20;

samp = 100; % number of PSD matrices
rep_A = 100; % number of iterations of PSD matrices
rep_B = 100; % number of iterations of perturbation term for each PSD matrices

dim = 50;

zeroprop = 0.2;

%% iteration of A

seed_A = 1;
ka = 87;
while ka < rep_A

    try

    % matrices before perturbation
    for i = 1:samp
        rng(10^6*seed_A + i)
        % random PSD matrices with close-to-zero eigenvalue(s)
        A(:,:,i) = PSD_defined_matrices(dim,round(zeroprop*dim),e1);
    end
    
    % means before perturbation
    mean_c1_bw = compute_Mean_cheap(A,eps,"BW"); % C
    mean_i1_bw = real(compute_Mean_inductive(A,eps,"BW")); % I
    mean_p1_bw = real(compute_Mean_projection2(A,eps1,eps2,seq,"BW")); % P

    ka = ka + 1; % number of valid iterations of B

    % distances
%     dist0(1,ka) = compute_F_norm(mean_i1_bw); % |I|
%     dist0(2,ka) = compute_F_norm(mean_p1_bw); % |P|
%     dist0(3,ka) = compute_F_norm(mean_c1_bw); % |C|
    dist0(4,ka) = compute_F_distance(mean_i1_bw, mean_p1_bw); % d(I,P)
    dist0(5,ka) = compute_F_distance(mean_c1_bw, mean_p1_bw); % d(C,P)
%     dist0(6,ka) = max(eig(mean_p1_bw)); % largest eigenvalue of P

    catch

    end

    seed_A = seed_A + 1; % track seed, can be >100

%% iteration of B
    dist = zeros(3);
    rob = zeros(1);

    seed_B = 1;
    kb = 0;
    while kb < rep_B
        
        try

        % matrices after perturbation
        for i = 1:samp
            rng(10^6*seed_A + 10^3*seed_B + i)
            % random perturbation terms
            b = randn(dim);
            Ax(:,:,i) = A(:,:,i) + e2*(b'+b)/2;
            % polar decomposition
            [R U V] = poldecomp(Ax(:,:,i));
            Ax_rie(:,:,i) = U;
            Ax_bw(:,:,i) = (Ax(:,:,i)+U)/2;
            Axx_rie(:,:,i) = Ax_rie(:,:,i) + e3*eye(dim);
            Axx_bw(:,:,i) = Ax_bw(:,:,i) + e3*eye(dim);
        end

        % means after perturbation
        mean_c2_bw = compute_Mean_cheap(Axx_bw,eps,"BW"); % C'
        mean_i2_bw = real(compute_Mean_inductive(Axx_bw,eps,"BW")); % I'
%         mean_p2_bw = real(compute_Mean_projection2(Axx_bw,eps1,eps2,seq,"BW")); % P'

        kb = kb + 1; % number of valid iterations of B

        % distances
%         dist(1,kb) = compute_F_distance(mean_i2_bw, mean_i1_bw); % d(I',I)
%         dist(2,kb) = compute_F_distance(mean_p2_bw, mean_p1_bw); % d(P',P)
%         dist(3,kb) = compute_F_distance(mean_c2_bw, mean_c1_bw); % d(C',C)
        dist(4,kb) = compute_F_distance(mean_i2_bw, mean_p1_bw); % d(I',P)
        dist(5,kb) = compute_F_distance(mean_c2_bw, mean_p1_bw); % d(C',P)

        catch

        end

        seed_B = seed_B + 1; % track seed, can be >100
        1000*(seed_A-1) + (seed_B-1) % show iteration progress

    end

    % output
    dist0_out = dist0; % 5*rep_A
    dist_out(:,:,ka) = dist; % 6*rep_B*rep_A

    save(sprintf('Results/2023/fdist0_dim%d_num%d.mat',dim,samp),'dist0_out')
    save(sprintf('Results/2023/fdist_dim%d_num%d.mat',dim,samp),'dist_out')
    
end


