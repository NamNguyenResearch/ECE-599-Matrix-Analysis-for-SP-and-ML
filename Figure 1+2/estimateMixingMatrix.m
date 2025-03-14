function [A_est,iter]=estimateMixingMatrix(SNR_dB,A,S)
    %====================== Time line of coding ===========================
    % Programmer: Nam Nguyen
    % Oregon State University, Corvallis, Oregon, United State
    % Date: May 25, 2022

    %====================== Usage =========================================

    % [A_est,iter]=estimateMixingMatrix(SNR_dB,A,S)

    %====================== Parameters ====================================
    % Outputs: 
    % A_est     :   estimated mixing matrix
    % iter      :   estimated iteration
    %----------------------------------------------------------------------
    % Inputs    :
    % SNR_dB    :   signal-to-noise ratio (SNR) in dB
    % A         :   mixing matrix
    % S         :   source signal
    %======================================================================
    
    
    %% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    global T;           % sequence length
    global N;           % number of sensors
    global K;           % number of sources
    global L_constant;  % fixed time frame
    global err_bound;   %tolerance between the currently estimated channel and the previously estimeated channel
    global iter_max;    % maximal number of iteration
    

    %% Calculate signal-to-noise ratio (SNR) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    buff=0;

    for t=1:T 
        buff=buff+mean((norm(A*S(:,t)))^2); 
    end

    buff=1/T*buff;

    SNR=10^(SNR_dB/10);

    variance=buff/SNR;
    
    
    %% Create the noise, v(t) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    mu=0;
    sigma=sqrt(variance);
    
    V=mu+sigma*randn(N,T);

    
    %% Create the received signal vector, x(t) %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    X=A*S+V;

    
    %% Estimate the matrix A by using the Prewhitened Alternating Projection Algorithm 
    [A_est,iter] = PAPA_algorithm(X,K,L_constant,err_bound,iter_max);
end