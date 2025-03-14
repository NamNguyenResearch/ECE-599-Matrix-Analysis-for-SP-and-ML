function [S]=createSourceSignal()
    %====================== Time line of coding ===========================
    % Programmer: Nam Nguyen
    % Oregon State University, Corvallis, Oregon, United State
    % Date: May 25, 2022

    %====================== Usage =========================================

    % [S]=createSourceSignal()

    %====================== Parameters ====================================
    % Outputs: 
    % S     :   source signal
    %======================================================================
    
    
    %% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    global L_low;   % the lowest value of local time frame
    global L_upp;   % the highest value of local time frame
    global T;       % sequence length
    global K;       % number of sources
    
    
    %% Generate the synthetic quasi-stationary singals %%%%%%%%%%%%%%%%%%%%
    S=zeros(K,T);
    S1=zeros(K,T);
    
    T_cur=1;

    while T_cur <= T
        for k=1:K
            L=round(unifrnd(L_low,L_upp));
            sigma_s=unifrnd(0,1);
        
            mu=0;                       % mean
            sigma=sqrt(sigma_s^2/2);    % standard deviation
            
            for temp=0:(L-1)
                temp1=T_cur+temp;
            
                s_R(temp1)=laprnd(1,1,mu,sigma);
                s_I(temp1)=laprnd(1,1,mu,sigma);
            
                S1(k,temp1)=s_R(temp1)+1i*s_I(temp1);
            end
        end
        T_cur=T_cur+L;
    end
    
    for k=1:K 
        for t=1:T 
           S(k,t)=S1(k,t); 
        end
    end
end