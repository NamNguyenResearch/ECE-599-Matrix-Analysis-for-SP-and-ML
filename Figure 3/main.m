%====================== Time line of coding ===============================
% Programmer: Nam Nguyen
% Oregon State University, Corvallis, Oregon, United State
% Date: May 25, 2022
%
% Thanks for the function "PAPA_mod.m" of authors 
% K. K. Lee, W.-K. Ma, X. Fu, T.-H. Chan, and C.-Y. Chi, "A Khatri-Rao
% Subspace Approach to Blind Identification of Mixtures of Quasi-Stationary
% Sources," Signal Processing, vol. 93, no. 12, pp. 3515-3527, Dec 2013
%=========================================================================

% Note: Please try to run this Matlab Code until getting the good figures (2-3 times)
% Thank you!

clear;
clc;


%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global L_low;       % the lowest value of local time frame
global L_upp;       % the highest value of local time frame
global T;           % sequence length
global N;           % number of sensors
global K;           % number of sources
global error_bound; % tolerance between the currently estimated channel and the previously estimeated channel
global iter_max;    % maximal number of iteration 
global SNR_dB;      % signal to noise ratio in dB

L_low=100;
L_upp=300;
T=79800;
N=5;
K=4; 
error_bound=10^(-6); 
iter_max=100;
SNR_dB=25;


%% Create the source signal vector, s(t) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
S=createSourceSignal();

%% Create the mixing matrix, A %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A=randn(N,K);
A=normalize(A,'norm',2);


%% Calculation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M=1:50:500;
MSE_dB=zeros(1,length(M));

for i=1:length(M)
    [MSE_dB(i)]=estimateMSE(M(i),A,S); 
end

    
%% Plot figure: MSE versus SNR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(3); 
plot(M,MSE_dB,'b--o','linewidth',1.5); 
grid on;
xlabel('Number of available frames, M','fontsize',16);
ylabel('MSE (dB)','fontsize',16);
title('The average MSE versus M','fontsize',16);