%Characterizations of the Delay Tail Distributions
%by Derya Malak
%Last revised: Nov 28, 2018

clear; clc; %close all


lineStyles = {'--','-','--','-.',':'};
markerStyles = {'none'};
colorStyles = {'k'};
pos = [950,200,500,380];

%Transmission time (\tau): # of frames being transmitted per successful frame

%Total delay (D): time from when a frame is first transmitted to when its ACK received

%# of channel states of a multistate Markov process S_t
%K = 1;      %Memoryless channel (There is only 1 state)
K = 2;      %There are Good G and Bad B states


r = 0.3;
rf = r;    %1/r represents the average error burst    %(or r = 0.3)    
rr = r;
if rf == 0
    K = 1;
end

scale = 1;
k = scale*5;      %round trip time (RTT)

eps_count = 20;
eps_set = linspace(0.001,0.5,eps_count) ; %the probability of block error (design requirement)
z = 1;     %This is the z-transform parameter

i = 20;
epsf = eps_set(i);
epsr = epsf;
eps = epsf;


T_count = 3;
T_set = scale*[15,10,8]; %timer=time to timeout (varies between k and 20)
%When timer (T) is increased, both throughput and delay are higher
%Throughput is upper bounded by (1-eps) as T->infty
%No such upper bound for delay  


T = T_set(1);


eps_Gf = 0;  %no error in good state
eps_Bf = 1;  %always error in bad state

eps_Gr = 0;  %no error in good state
eps_Br = 1;  %always error in bad state


%Coded ARQ Model
M = 2;      %total number of transmitted packets
N = 2;      %required number of successfully received packets   


NACK = 1; %1 with NACK and 0 without NACK

symbolic = 0; %1 for memoryless channel results

%Generate the probability transition matrices
[P00,P01,P10,P11,P0x,Px0,P1x,Px1,P_kron,pi_kron,pi_I_kron] = TransitionProbabilityMatrices(eps_Gf,eps_Bf,epsf,eps_Gr,eps_Br,epsr,rf,rr,NACK,K);

Ztransform = 1; %Take inverse Z transform to compute the PMFs and CCDFs of
%the delay (computationally costly)
%Ztransform = 0; %Use numerical inversion to approximate the PMFs and CCDFs of the delay (less complex)

%% RUN THE FOLLOWING SECTIONS INDIVIDUALLY TO COMPUTE THE DISTRIBUTION OF THE DELAY
% MEMORYLESS CHANNEL (I.E. WHEN K=1), WE COMPUTE THE DISTRIBUTION BY TAKING THE INVERSE Z TRANSFORM
% GILBERT ELLIOTT CHANNEL (I.E. WHEN K=2),  WE COMPUTE THE DISTRIBUTION BY
% TAKING THE INVERSE OF THE CHARACTERISTIC FUNCTION NUMERICALLY

%% 1-Uncoded ARQ

%Use matrix-generating function of delay to plot the delay PMF (uncoded ARQ)
[phiD_ARQ, PMF_delay_ARQ, CCDF_delay_ARQ] = UncodedARQtailDelay(k,K,T,P_kron,pi_I_kron,P00,P01,P10,P11,P0x,Px0,P1x,Px1,epsf,Ztransform);

%% 2-Uncoded HARQ

%Use matrix-generating function of delay to plot the delay PMF (uncoded HARQ)
[phiD_HARQ, PMF_delay_HARQ, CCDF_delay_HARQ] = HARQtailDelay(k,K,T,P_kron,pi_I_kron,P00,P01,P10,P11,P0x,Px0,P1x,Px1,r,NACK,epsf,eps,Ztransform);

%% 3-CF ARQ

%Use matrix-generating function of delay to plot the delay PMF (CF-ARQ)
[phiD_CF_ARQ, PMF_delay_CF_ARQ, CCDF_delay_CF_ARQ] = CF_ARQtailDelay(k,K,T,P_kron,pi_I_kron,P00,P01,P10,P11,P0x,Px0,P1x,Px1,epsf,M,Ztransform);


%% 4-Coded ARQ

%Use matrix-generating function of delay to plot the delay PMF (C-ARQ)
[phiD_C_ARQ, PMF_delay_C_ARQ, CCDF_delay_C_ARQ] = C_ARQtailDelay(k,K,T,P_kron,pi_I_kron,P00,P01,P10,P11,P0x,Px0,P1x,Px1,epsf,M,Ztransform);

