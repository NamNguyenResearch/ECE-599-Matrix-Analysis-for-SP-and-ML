function [A_est,iter]=PAPA_algorithm(X,K,L,error_bound,iter_max)
    %====================== Time line of coding ===========================
    % Programmer: Nam Nguyen
    % Oregon State University, Corvallis, Oregon, United State
    % Date: May 25, 2022

    %====================== Usage =========================================

    % [A_est,iter]=PAPA_algorithm(X,K,L,error_bound,iter_max)

    %====================== Parameters ====================================
    % Outputs: 
    % A_est         :   estimated mixing matrix
    % iter          :   estimated iteration 
    %----------------------------------------------------------------------
    % Inputs:
    % X             :   received singal
    % K             :   number of sources
    % L             :   frame length
    % error_bound   :   tolerance between the currently estimated channel and the previously estimeated channel
    % iter_max      :   maximal number of iteration 
    %======================================================================
    
    
    %% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    global T;           % sequence length
    global N;           % number of sensors

    
    %% Calculate local covariance %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    M=floor(T/L);       % Number of frames
    M_overlap=2*M-1;    % number of frames after overlapping
    
    for temp=1:M-1
        R(:,:,temp)=X(:,(temp-1)*L+1:temp*L)*X(:,(temp-1)*L+1:temp*L)'/L;
        R(:,:,M+temp)=X(:,(temp-1)*L+1+(L/2):temp*L+(L/2))*X(:,(temp-1)*L+1+(L/2):temp*L+(L/2))'/L;
    end
    
    R(:,:,M)=X(:,(M-1)*L+1:M*L)*X(:,(M-1)*L+1:M*L)'/L;


    %% Remove noise %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    sigma_est=zeros(M_overlap,1);
    
    % Estimate sigma %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for temp=1:M_overlap
        [eigvec,eigval]=eig(R(:,:,temp));
        if min(diag(eigval))<sigma_est
            sigma_est(temp)=min(diag(eigval));
        end
    end
    
    % Subtract noise %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for temp=1:M_overlap
        R_noisefree(:,:,temp)=R(:,:,temp)-min(sigma_est)*eye(N);
    end

    
    %% Pre-whitening process %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    R=R_noisefree;

    [N,N,M]=size(R);
    
    % Time-averaged global covariance %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    R_bar_vec=reshape(R,N^2,M);
    R_bar_vec=(1/M)*sum(R_bar_vec,2); 
    R_bar=reshape(R_bar_vec,N,N);
    
    % Square-root factorization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [U,S]=svd(0.5*(R_bar+R_bar'));
    B=U(:,1:K)*sqrt(S(1:K,1:K));
   
    PinvB=pinv(B);                              % pseudo inverse of B 
    B_pre=B;
    
    R_vec=reshape(R,N^2,M);
    R_whiten=kron(conj(PinvB),PinvB)*R_vec;     % prewhitening
   

    %% Khatri-Rao subspace extraction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tic
    
    z=R_whiten(:);
    Y=reshape(z,K^2,M);

    iter=0;
    
    C=Y*Y';
    [U,S]=eig(C+C');

    Us=U(:,end-K+1:end);

    for k=1:K
        % Random initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        h=Us*randn(K,1);   
    
        % Alternating Projection %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [g,g_Kronecker,iter_t]=AP_algorithm(h,iter_max,error_bound,Us,K);  
        iter=iter+iter_t;
        G(:,k)=g;

        if k < K
            % Orthogonal projection %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            U_temp=g_Kronecker'*Us;
            Us=Us-g_Kronecker*U_temp;
        end
    end

    time_PAPA=toc;
    
    iter=iter/K; 
    A_est=B_pre*G;
end