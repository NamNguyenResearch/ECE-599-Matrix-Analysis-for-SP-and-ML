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
global L_constant;  % fixed time frame
global err_bound;   % tolerance between the currently estimated channel and the previously estimeated channel
global iter_max;    % maximal number of iteration 

L_low=100;
L_upp=300;
T=79800;
N=5;
K=4; 
L_constant=200;
err_bound=10^(-6); 
iter_max=100;


%% Create the source signal vector, s(t) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S=createSourceSignal();


%% Create the mixing matrix, A %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A=randn(N,K);
A=normalize(A,'norm',2);


%% Calculation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SNR_dB=-10:5:40;
MSE_dB=zeros(1,length(SNR_dB));

for i=1:length(SNR_dB)
    [MSE_dB(i)]=estimateMSE(SNR_dB(i),A,S); 
end


%% Plot figure: Segment of synthetic quasi-stationary source sginals %%%%%%
figure(1)
subplot(4,1,1)
plot(1:7000,real(S(1,1:7000)),'b-','linewidth',1.5); 
    
subplot(4,1,2)
plot(1:7000,real(S(2,1:7000)),'b-','linewidth',1.5);
    
subplot(4,1,3)
plot(1:7000,real(S(3,1:7000)),'b-','linewidth',1.5);
    
subplot(4,1,4)
plot(1:7000,real(S(4,1:7000)),'b-','linewidth',1.5);
    
    
%% Plot figure: MSE versus SNR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(2); 
plot(SNR_dB,MSE_dB,'b--o','linewidth',1.5); 
grid on;
xlabel('SNR (dB)','fontsize',16);
ylabel('MSE (dB)','fontsize',16);
title('The average MSE versus SNR','fontsize',16);