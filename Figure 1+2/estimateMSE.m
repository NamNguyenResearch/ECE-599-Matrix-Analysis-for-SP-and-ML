function [MSE_dB]=estimateMSE(SNR_dB,A,S)
    %====================== Time line of coding ===========================
    % Programmer: Nam Nguyen
    % Oregon State University, Corvallis, Oregon, United State
    % Date: May 25, 2022

    %====================== Usage =========================================

    % [MSE_dB]=estimateMSE(SNR_dB,A,S)

    %====================== Parameters ====================================
    % Outputs: 
    % MSE_dB    :   average mean square error (MSE)
    %----------------------------------------------------------------------
    % Inputs    :
    % SNR_dB    :   signal-to-noise ratio (SNR) in dB
    % A         :   mixing matrix
    % S         :   source signal
    %======================================================================
    
    
    %% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    global K;    % number of sources
    
    
    %% Estimate mixing matrix and iteration %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [A_est, iter]=estimateMixingMatrix(SNR_dB,A,S);
    

    %% Calculate Average mean square error (MSE) %%%%%%%%%%%%%%%%%%%%%%%%%%
    for k=1:K   
        a=normalize(A(:,k),'norm',2);
   
        for j=1:K  
            b=normalize(A_est(:,j),'norm',2);
             
            Sub1(k,j)=(norm(a-(1)*b))^2; 
            Sub2(k,j)=(norm(a-(-1)*b))^2;
            
            Sub1(k,j+K)=Sub2(k,j);
        end
    end

    temp=1;  
    

    for i=1:(2*K) 
        for j=1:(2*K) 
            for m=1:(2*K) 
                for n=1:(2*K) 
                    Sum(temp)=Sub1(1,i)+Sub1(2,j)+Sub1(3,m)+Sub1(4,n);  
                    temp=temp+1;
                end
            end
        end
    end

    MSE=min(Sum/K);
    MSE_dB=10*log10(MSE);
end