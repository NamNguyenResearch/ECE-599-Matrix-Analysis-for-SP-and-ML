function [y] = laprnd(m, n, mu, sigma)
    %====================== Time line of coding ===========================
    % Programmer: Nam Nguyen
    % Oregon State University, Corvallis, Oregon, United State
    % Date: May 25, 2022

    %====================== Usage =========================================

    % [y] = laprnd(m, n, mu, sigma)

    %====================== Parameters ====================================
    % Outputs: 
    % y     :   laplacian random number drawn from Laplacian distribution
    %----------------------------------------------------------------------
    % Inputs:
    % m     :   number of rows
    % n     :   number of columns
    % mu    :   mean
    % sigma :   standard deviation
    %======================================================================
    
    
    %% Calculation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    u=rand(m,n)-0.5;
    b=sigma/sqrt(2);
    y=mu-b*sign(u).*log(1-2*abs(u));
end
