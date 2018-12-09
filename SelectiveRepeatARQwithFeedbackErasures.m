% by Derya Malak
% Last revised: Feb 27, 2018 Tuesday
% ARQ model with retransmissions and feedback erasures (for uncoded and coded schemes)
% and/or HARQ with soft combining

clear; clc; %close all; 

lineStyles = {'--','-','--','-.',':'};
markerStyles = {'none'};
colorStyles = {'k'};
pos = [950,200,500,380];

%Transmission time (\tau): # of frames being transmitted per successful frame

%Total delay (D): time from when a frame is first transmitted to when its ACK received

%# of channel states of a multistate Markov process S_t
K = 2;      %There are Good G and Bad B states

r = 0.3; 
%r = 0.1;
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

T_count = 3;
T_set = scale*[15,10,8]; %timer=time to timeout (varies between k and 20)
%When timer (T) is increased, both throughput and delay are higher
%Throughput is upper bounded by (1-eps) as T->infty
%No such upper bound for delay     


eps_Gf = 0;  %no error in good state
eps_Bf = 1;  %always error in bad state

eps_Gr = 0;  %no error in good state
eps_Br = 1;  %always error in bad state


%Coded ARQ Model
M = 2;      %total number of transmitted packets
N = 2;      %required number of successfully received packets   


NACK = 1; %1 with NACK and 0 without NACK

symbolic = 0; %1 for memoryless channel results

%% RUN THE FOLLOWING SECTIONS INDIVIDUALLY TO CHARACTERIZE THE PERFORMANCE OF THE DELAY AND THE THROUGHPUT
% MEMORYLESS CHANNEL (I.E. WHEN K=1)
% GILBERT ELLIOTT CHANNEL (I.E. WHEN K=2)
% THE FUNCTION "TransitionProbabilityMatrices" COMPUTES THE TRANSITION
% PROBABILITY MATRICES OF THE CHANNEL

% THE FOLLOWING FUNCTIONS RETURN THE DELAY PGF, AVERAGE DELAY, VARIANCE OF
% THE DELAY, AND THROUGHPUT FOR THE RESPECTIVE SCHEMES

% THE FUNCTION "NoCodingPhiD" FOR UNCODED ARQ
% THE FUNCTION "NoCodingHARQPhiD" FOR UNCODED HARQ
% THE FUNCTION "CumulativeFeedbackPhiD" FOR CUMULATIVE FEEDBACK ARQ
% THE FUNCTION "CodedPhiD" FOR CODED ARQ

%% 1-Uncoded ARQ scheme with no HARQ combining

meanDelay1 = zeros(T_count,eps_count);
Throughput1 = zeros(T_count,eps_count);

varDelay1 = zeros(T_count,eps_count);

%DELAY vs block-error rate \epsilon
figure
clear str;      str = cell(1,T_count);
for t = 1:T_count
    T = T_set(t);
    str{t} = ['Uncoded ARQ, T=' num2str(T)]; 
    
    for i = 1:eps_count
        epsf = eps_set(i);
        epsr = epsf;

        [P00,P01,P10,P11,P0x,Px0,P1x,Px1,P_kron,pi_kron,pi_I_kron] = TransitionProbabilityMatrices(eps_Gf,eps_Bf,epsf,eps_Gr,eps_Br,epsr,rf,rr,NACK,K);
        [PhiD1, phiD1, meanDelay1(t,i), varDelay1(t,i), Throughput1(t,i)] = NoCodingPhiD(z,k,K,T,P_kron,pi_I_kron,P00,P01,P10,P11,Px0,Px1,symbolic);  
    end
    plot(eps_set,meanDelay1(t,:),'k','linewidth',2,...
        'color',colorStyles{1+rem(t,numel(colorStyles))},'linestyle',lineStyles{1+rem(t,numel(lineStyles))},...
        'marker',markerStyles{1+rem(t,numel(markerStyles))},'markersize',12); hold on;
end
xlab = 'Erasure rate, $\epsilon$'; 
ylab = 'Average delay, $\bar{D}$';
box on;     set(gca,'FontSize',20) 
xlhand = get(gca,'xlabel'); xlabel(xlab,'Interpreter','latex'); set(xlhand,'fontsize',20) 
ylhand = get(gca,'ylabel'); ylabel(ylab,'Interpreter','latex'); set(ylhand,'fontsize',20)
legend(str,'FontSize',20,'location','NorthEast');
xaxis = 0.3;    yaxis = 5;
text(xaxis,yaxis,['r=' num2str(r) ', k=' num2str(k)],'fontsize',20);



%THROUGHPUT vs block-error rate \epsilon
figure
for t = 1:T_count
    T = T_set(t);
     
    plot(eps_set,Throughput1(t,:),'k','linewidth',2,...
        'color',colorStyles{1+rem(t,numel(colorStyles))},'linestyle',lineStyles{1+rem(t,numel(lineStyles))},...
        'marker',markerStyles{1+rem(t,numel(markerStyles))},'markersize',12); hold on;
end
xlab = 'Erasure rate, $\epsilon$'; 
ylab = 'Throughput, $\eta$';
box on;     set(gca,'FontSize',20) 
xlhand = get(gca,'xlabel'); xlabel(xlab,'Interpreter','latex'); set(xlhand,'fontsize',20) 
ylhand = get(gca,'ylabel'); ylabel(ylab,'Interpreter','latex'); set(ylhand,'fontsize',20)
legend(str,'FontSize',20,'location','NorthEast');
xaxis = 0.3;    yaxis = 0.35;
text(xaxis,yaxis,['r=' num2str(r) ', k=' num2str(k)],'fontsize',20);

%{
%DELAY SECOND MOMENT vs block-error rate \epsilon
figure
for t = 1:T_count
    T = T_set(t);
    
    plot(eps_set,varDelay1(t,:)+meanDelay1(t,:).^2,'k','linewidth',2,...
        'color',colorStyles{1+rem(t,numel(colorStyles))},'linestyle',lineStyles{1+rem(t,numel(lineStyles))},...
        'marker',markerStyles{1+rem(t,numel(markerStyles))},'markersize',12); hold on;
end
xlab = 'Erasure rate, $\epsilon$'; 
ylab = 'Second moment of delay, $E[D^2]$';
box on;     set(gca,'FontSize',20) 
xlhand = get(gca,'xlabel'); xlabel(xlab,'Interpreter','latex'); set(xlhand,'fontsize',20) 
ylhand = get(gca,'ylabel'); ylabel(ylab,'Interpreter','latex'); set(ylhand,'fontsize',20)
legend(str,'FontSize',20,'location','NorthEast');
xaxis = 0.3;    yaxis = 5;
text(xaxis,yaxis,['r=' num2str(r) ', k=' num2str(k)],'fontsize',20);
%}

%DELAY VARIABILITY vs block-error rate \epsilon
figure
for t = 1:T_count
    T = T_set(t);
    
    plot(eps_set,varDelay1(t,:),'k','linewidth',2,...
        'color',colorStyles{1+rem(t,numel(colorStyles))},'linestyle',lineStyles{1+rem(t,numel(lineStyles))},...
        'marker',markerStyles{1+rem(t,numel(markerStyles))},'markersize',12); hold on;
end
xlab = 'Erasure rate, $\epsilon$'; 
ylab = 'Variability of delay, $\sigma^2_D$';
box on;     set(gca,'FontSize',20) 
xlhand = get(gca,'xlabel'); xlabel(xlab,'Interpreter','latex'); set(xlhand,'fontsize',20) 
ylhand = get(gca,'ylabel'); ylabel(ylab,'Interpreter','latex'); set(ylhand,'fontsize',20)
legend(str,'FontSize',20,'location','NorthEast');
xaxis = 0.3;    yaxis = 5;
text(xaxis,yaxis,['r=' num2str(r) ', k=' num2str(k)],'fontsize',20);


%Guaranteeable DELAY vs block-error rate \epsilon
figure
for t = 1:T_count
    T = T_set(t);
    
    plot(eps_set,meanDelay1(t,:)+2*sqrt(varDelay1(t,:)),'k','linewidth',2,...
        'color',colorStyles{1+rem(t,numel(colorStyles))},'linestyle',lineStyles{1+rem(t,numel(lineStyles))},...
        'marker',markerStyles{1+rem(t,numel(markerStyles))},'markersize',12); hold on;
end
xlab = 'Erasure rate, $\epsilon$'; 
ylab = 'Guaranteeable delay, $\hat{D}$';
box on;     set(gca,'FontSize',20) 
xlhand = get(gca,'xlabel'); xlabel(xlab,'Interpreter','latex'); set(xlhand,'fontsize',20) 
ylhand = get(gca,'ylabel'); ylabel(ylab,'Interpreter','latex'); set(ylhand,'fontsize',20)
legend(str,'FontSize',20,'location','NorthEast');
xaxis = 0.3;    yaxis = 5;
text(xaxis,yaxis,['r=' num2str(r) ', k=' num2str(k)],'fontsize',20);


%DELAY VARIABILITY vs RTT k      
kcount = 20;
kset = linspace(.1,10,kcount); %values for round trip time (RTT)


count = 4;
meanDelay1 = zeros(count,kcount);
Throughput1 = zeros(count,kcount);

varDelay1 = zeros(count,kcount);

figure
clear str;      str = cell(1,count);
count = 0 ;
for eps = [0, 0.1, 0.25, 0.5]     
    count = count+1;
    str{count} = ['Uncoded ARQ, \epsilon=' num2str(eps)]; 
    
    for kval = 1:kcount
        
        k = kset(kval);
        T = kval*kset(kcount);       

        [P00,P01,P10,P11,P0x,Px0,P1x,Px1,P_kron,pi_kron,pi_I_kron] = TransitionProbabilityMatrices(eps_Gf,eps_Bf,epsf,eps_Gr,eps_Br,epsr,rf,rr,NACK,K);
        [PhiD1, phiD1, meanDelay1(count,kval), varDelay1(count,kval), Throughput1(count,kval)] = NoCodingPhiD(z,k,K,T,P_kron,pi_I_kron,P00,P01,P10,P11,Px0,Px1,symbolic);  
    end
    plot(kset/T,varDelay1(count,:),'k','linewidth',2,...
        'color',colorStyles{1+rem(t,numel(colorStyles))},'linestyle',lineStyles{1+rem(t,numel(lineStyles))},...
        'marker',markerStyles{1+rem(t,numel(markerStyles))},'markersize',12); hold on;
end
xlab = '$k/T$'; 
ylab = 'Variability of delay, $\sigma^2_D$';
box on;     set(gca,'FontSize',20) 
xlhand = get(gca,'xlabel'); xlabel(xlab,'Interpreter','latex'); set(xlhand,'fontsize',20) 
ylhand = get(gca,'ylabel'); ylabel(ylab,'Interpreter','latex'); set(ylhand,'fontsize',20)
legend(str,'FontSize',20,'location','NorthEast');
xaxis = 0.3;    yaxis = 5;
text(xaxis,yaxis,['r=' num2str(r) ', T=2k'],'fontsize',20);


k = scale*5;


%{
pi_I, pi_B and pi_C are the probability vectors of states I, B, and C
%b*inv(A) = b/A         %inv(A)*b = A\b


A1 = eye(2*K,2*K)-P_kron^(k-1)*(P10+P11*P_kron^d);
B1 = pi_I_kron*P_kron^(k-1)*(P10+P11*P_kron^d);
pi_B_kron = B1/A1;      


A2 = eye(2*K,2*K)-Px1^T;
B2 = (pi_I_kron+pi_B_kron)*P_kron^(k-1)*P01*Px1^d;
pi_C_kron = B2/A2;      


pi_kron_est = pi_I_kron + pi_B_kron + pi_C_kron;


pi_kron - pi_kron_est
sum(pi_kron)
sum(pi_kron_est)
%}


%% 2-Uncoded ARQ scheme and HARQ with soft combining
%%{
meanDelay2 = zeros(T_count,eps_count);
Throughput2 = zeros(T_count,eps_count);


%DELAY vs block-error rate \epsilon
figure
clear str;      str = cell(1,T_count);
for t = 1:T_count
    T = T_set(t);
    str{t} = ['HARQ, T=' num2str(T)]; 
    
    for i = 1:eps_count
        eps = eps_set(i);
        %[P00,P01,P10,P11,P0x,Px0,P1x,Px1,P_kron,pi_kron,pi_I_kron] = TransitionProbabilityMatrices(eps_Gf,eps_Bf,epsf,eps_Gr,eps_Br,epsr,rf,rr,NACK,K);
        [PhiD2, phiD2, meanDelay2(t,i), Throughput2(t,i)] = NoCodingHARQPhiD(z,k,K,T,eps,r,NACK,symbolic);                                                             
    end
    plot(eps_set,meanDelay2(t,:),'k','linewidth',2,...
        'color',colorStyles{1+rem(t,numel(colorStyles))},'linestyle',lineStyles{1+rem(t,numel(lineStyles))},...
        'marker',markerStyles{1+rem(t,numel(markerStyles))},'markersize',12); hold on;
end
xlab = 'Erasure rate, $\epsilon$'; 
ylab = 'Average delay, $\bar{D}$';
box on;     set(gca,'FontSize',20) 
xlhand = get(gca,'xlabel'); xlabel(xlab,'Interpreter','latex'); set(xlhand,'fontsize',20) 
ylhand = get(gca,'ylabel'); ylabel(ylab,'Interpreter','latex'); set(ylhand,'fontsize',20)
legend(str,'FontSize',20,'location','NorthEast');
xaxis = 0.3;    yaxis = 5;
text(xaxis,yaxis,['r=' num2str(r) ', k=' num2str(k)],'fontsize',20);


%THROUGHPUT vs block-error rate \epsilon
figure
for t = 1:T_count
    T = T_set(t);
     
    plot(eps_set,Throughput2(t,:),'k','linewidth',2,...
        'color',colorStyles{1+rem(t,numel(colorStyles))},'linestyle',lineStyles{1+rem(t,numel(lineStyles))},...
        'marker',markerStyles{1+rem(t,numel(markerStyles))},'markersize',12); hold on;
end
xlab = 'Erasure rate, $\epsilon$'; 
ylab = 'Throughput, $\eta$';
box on;     set(gca,'FontSize',20) 
xlhand = get(gca,'xlabel'); xlabel(xlab,'Interpreter','latex'); set(xlhand,'fontsize',20) 
ylhand = get(gca,'ylabel'); ylabel(ylab,'Interpreter','latex'); set(ylhand,'fontsize',20)
legend(str,'FontSize',20,'location','NorthEast');
xaxis = 0.3;    yaxis = 0.35;
text(xaxis,yaxis,['r=' num2str(r) ', k=' num2str(k)],'fontsize',20);
%}


%% 3-Cumulative feedback-based ARQ (CF ARQ) scheme with no HARQ combining 

meanDelay3 = zeros(T_count,eps_count);
Throughput3 = zeros(T_count,eps_count);


%DELAY (average per packet) vs block-error rate \epsilon
figure
clear str;      str = cell(1,T_count);
for t = 1:T_count
    T = T_set(t);
    str{t} = ['CF ARQ, T=' num2str(T)]; 
    
    for i = 1:eps_count
        epsf = eps_set(i);
        epsf = sqrt(epsf^4+2*epsf^3*(1-epsf));        
        epsr = epsf;
        
        if rf == 0
           eps_Bf = epsf;
        end
        
        if rr == 0
           eps_Br = epsr;
        end
          
        [P00,P01,P10,P11,P0x,Px0,P1x,Px1,P_kron,pi_kron,pi_I_kron] = TransitionProbabilityMatrices(eps_Gf,eps_Bf,epsf,eps_Gr,eps_Br,epsr,rf,rr,NACK,K);
        [PhiD3, phiD3, meanDelay3(t,i), Throughput3(t,i)] = CumulativeFeedbackPhiDv2(z,k,K,T,P_kron,pi_I_kron,P00,P01,P10,P11,Px0,Px1,M,N,symbolic);

    end
    plot(eps_set,meanDelay3(t,:),'k','linewidth',2,...
        'color',colorStyles{1+rem(t,numel(colorStyles))},'linestyle',lineStyles{1+rem(t,numel(lineStyles))},...
        'marker',markerStyles{1+rem(t,numel(markerStyles))},'markersize',12); hold on;
end
xlab = 'Erasure rate, $\epsilon$'; 
ylab = 'Average delay, $\bar{D}$';
box on;     set(gca,'FontSize',20) 
xlhand = get(gca,'xlabel'); xlabel(xlab,'Interpreter','latex'); set(xlhand,'fontsize',20) 
ylhand = get(gca,'ylabel'); ylabel(ylab,'Interpreter','latex'); set(ylhand,'fontsize',20)
legend(str,'FontSize',20,'location','NorthEast');
xaxis = 0.3;    yaxis = 2;
text(xaxis,yaxis,['M=' num2str(M) ', N=' num2str(N) ', r=' num2str(r) ', k=' num2str(k)],'fontsize',20);

 
%THROUGHPUT vs block-error rate \epsilon
figure
for t = 1:T_count
    T = T_set(t);
      
     plot(eps_set,Throughput3(t,:),'k','linewidth',2,...
        'color',colorStyles{1+rem(t,numel(colorStyles))},'linestyle',lineStyles{1+rem(t,numel(lineStyles))},...
        'marker',markerStyles{1+rem(t,numel(markerStyles))},'markersize',12); hold on;
end
xlab = 'Erasure rate, $\epsilon$'; 
ylab = 'Throughput, $\eta$';
box on;     set(gca,'FontSize',20) 
xlhand = get(gca,'xlabel'); xlabel(xlab,'Interpreter','latex'); set(xlhand,'fontsize',20) 
ylhand = get(gca,'ylabel'); ylabel(ylab,'Interpreter','latex'); set(ylhand,'fontsize',20)
legend(str,'FontSize',20,'location','NorthEast');
xaxis = 0.3;    yaxis = 2;
text(xaxis,yaxis,['M=' num2str(M) ', N=' num2str(N) ', r=' num2str(r) ', k=' num2str(k)],'fontsize',20);


%% 4-Coded ARQ scheme with no HARQ combining 

meanDelay4 = zeros(T_count,eps_count);
Throughput4 = zeros(T_count,eps_count);
varDelay4 = zeros(T_count,eps_count);

%DELAY (average per packet) vs block-error rate \epsilon
figure
clear str;      str = cell(1,T_count);
for t = 1:T_count
    T = T_set(t);
    str{t} = ['Coded ARQ, T=' num2str(T)]; 
    
    for i = 1:eps_count
        epsf = eps_set(i);
        epsf = sqrt(epsf^4+2*epsf^3*(1-epsf));
        
        if rf == 0
           eps_Bf = epsf;
        end
                       
        epsr = epsf;
        
        if rr == 0
           eps_Br = epsr;
        end
           
        [P00,P01,P10,P11,P0x,Px0,P1x,Px1,P_kron,pi_kron,pi_I_kron] = TransitionProbabilityMatrices(eps_Gf,eps_Bf,epsf,eps_Gr,eps_Br,epsr,rf,rr,NACK,K);
       %[PhiD3, phiD3, meanDelay3(t,i), Throughput3(t,i)] = PipelinePhiD(z,k,K,T,P_kron,pi_I_kron,P00,P01,P10,P11,Px0,Px1,M,N);
       %[PhiD3, phiD3, meanDelay3(t,i), Throughput3(t,i)] = CumulativeFeedbackPhiD(z,k,K,T,P_kron,pi_I_kron,P00,P01,P10,P11,Px0,Px1,M,N,symbolic);
       [PhiD4, phiD4, meanDelay4(t,i), varDelay4(t,i), Throughput4(t,i)] = CodedPhiDv2(z,k,K,T,P_kron,pi_I_kron,P00,P01,P10,P11,Px0,Px1,M,N,symbolic);

    end
    plot(eps_set,meanDelay4(t,:),'k','linewidth',2,...
        'color',colorStyles{1+rem(t,numel(colorStyles))},'linestyle',lineStyles{1+rem(t,numel(lineStyles))},...
        'marker',markerStyles{1+rem(t,numel(markerStyles))},'markersize',12); hold on;
end
xlab = 'Erasure rate, $\epsilon$'; 
ylab = 'Average delay, $\bar{D}$';
box on;     set(gca,'FontSize',20) 
xlhand = get(gca,'xlabel'); xlabel(xlab,'Interpreter','latex'); set(xlhand,'fontsize',20) 
ylhand = get(gca,'ylabel'); ylabel(ylab,'Interpreter','latex'); set(ylhand,'fontsize',20)
legend(str,'FontSize',20,'location','NorthEast');
xaxis = 0.3;    yaxis = 2;
text(xaxis,yaxis,['M=' num2str(M) ', N=' num2str(N) ', r=' num2str(r) ', k=' num2str(k)],'fontsize',20);

 
%THROUGHPUT vs block-error rate \epsilon
figure
for t = 1:T_count
    T = T_set(t);
    
    plot(eps_set,M*Throughput4(t,:),'k','linewidth',2,...
        'color',colorStyles{1+rem(t,numel(colorStyles))},'linestyle',lineStyles{1+rem(t,numel(lineStyles))},...
        'marker',markerStyles{1+rem(t,numel(markerStyles))},'markersize',12); hold on;
end
xlab = 'Erasure rate, $\epsilon$'; 
ylab = 'Throughput, $\eta$';
box on;     set(gca,'FontSize',20) 
xlhand = get(gca,'xlabel'); xlabel(xlab,'Interpreter','latex'); set(xlhand,'fontsize',20) 
ylhand = get(gca,'ylabel'); ylabel(ylab,'Interpreter','latex'); set(ylhand,'fontsize',20)
legend(str,'FontSize',20,'location','NorthEast');
xaxis = 0.3;    yaxis = 2;
text(xaxis,yaxis,['M=' num2str(M) ', N=' num2str(N) ', r=' num2str(r) ', k=' num2str(k)],'fontsize',20);


%DELAY VARIABILITY vs block-error rate \epsilon
figure
for t = 1:T_count
    T = T_set(t);
    
    plot(eps_set,varDelay4(t,:),'k','linewidth',2,...
        'color',colorStyles{1+rem(t,numel(colorStyles))},'linestyle',lineStyles{1+rem(t,numel(lineStyles))},...
        'marker',markerStyles{1+rem(t,numel(markerStyles))},'markersize',12); hold on;
end
xlab = 'Erasure rate, $\epsilon$'; 
ylab = 'Variability of delay, $\sigma^2_D$';
box on;     set(gca,'FontSize',20) 
xlhand = get(gca,'xlabel'); xlabel(xlab,'Interpreter','latex'); set(xlhand,'fontsize',20) 
ylhand = get(gca,'ylabel'); ylabel(ylab,'Interpreter','latex'); set(ylhand,'fontsize',20)
legend(str,'FontSize',20,'location','NorthEast');
xaxis = 0.3;    yaxis = 5;
text(xaxis,yaxis,['r=' num2str(r) ', k=' num2str(k)],'fontsize',20);


%Guaranteeable DELAY vs block-error rate \epsilon
figure
for t = 1:T_count
    T = T_set(t);
    
    plot(eps_set,meanDelay4(t,:)+3*sqrt(varDelay4(t,:)),'k','linewidth',2,...
        'color',colorStyles{1+rem(t,numel(colorStyles))},'linestyle',lineStyles{1+rem(t,numel(lineStyles))},...
        'marker',markerStyles{1+rem(t,numel(markerStyles))},'markersize',12); hold on;
end
xlab = 'Erasure rate, $\epsilon$'; 
ylab = 'Guaranteeable delay, $\hat{D}$';
box on;     set(gca,'FontSize',20) 
xlhand = get(gca,'xlabel'); xlabel(xlab,'Interpreter','latex'); set(xlhand,'fontsize',20) 
ylhand = get(gca,'ylabel'); ylabel(ylab,'Interpreter','latex'); set(ylhand,'fontsize',20)
legend(str,'FontSize',20,'location','NorthEast');
xaxis = 0.3;    yaxis = 5;
text(xaxis,yaxis,['r=' num2str(r) ', k=' num2str(k)],'fontsize',20);
%legend([5,6,7,8],'Coded ARQ, T=15','Coded ARQ, T=8','Uncoded ARQ, T=15','Uncoded ARQ, T=8')



%DELAY VARIABILITY vs RTT k      
kcount = 20;
kset = linspace(.1,10,kcount); %values for round trip time (RTT)

count = 4;
meanDelay4 = zeros(count,kcount);
Throughput4 = zeros(count,kcount);

varDelay4 = zeros(count,kcount);

figure
clear str;      str = cell(1,count);
count = 0 ;
for eps = [0, 0.1, 0.25, 0.5]     
    count = count+1;
    str{count} = ['Coded ARQ, \epsilon=' num2str(eps)]; 
    
    for kval = 1:kcount
        
        k = kset(kval);
        T = kval*kset(kcount);

        [P00,P01,P10,P11,P0x,Px0,P1x,Px1,P_kron,pi_kron,pi_I_kron] = TransitionProbabilityMatrices(eps_Gf,eps_Bf,epsf,eps_Gr,eps_Br,epsr,rf,rr,NACK,K);
        [PhiD4, phiD4, meanDelay4(count,kval), varDelay4(count,kval), Throughput4(count,kval)] = CodedPhiD(z,k,K,T,P_kron,pi_I_kron,P00,P01,P10,P11,Px0,Px1,M,N,symbolic);                                                                                                                
    end
    plot(kset/T,varDelay4(count,:),'k','linewidth',2,...
        'color',colorStyles{1+rem(t,numel(colorStyles))},'linestyle',lineStyles{1+rem(t,numel(lineStyles))},...
        'marker',markerStyles{1+rem(t,numel(markerStyles))},'markersize',12); hold on;
end
xlab = '$k/T$'; 
ylab = 'Variability of delay, $\sigma^2_D$';
box on;     set(gca,'FontSize',20) 
xlhand = get(gca,'xlabel'); xlabel(xlab,'Interpreter','latex'); set(xlhand,'fontsize',20) 
ylhand = get(gca,'ylabel'); ylabel(ylab,'Interpreter','latex'); set(ylhand,'fontsize',20)
legend(str,'FontSize',20,'location','NorthEast');
xaxis = 0.3;    yaxis = 5;
text(xaxis,yaxis,['r=' num2str(r) ', T=2k'],'fontsize',20);


k = scale*5;


%% OTHER RELATED PLOTS

%Throughput and delay versus timeout T

markerStyles = {'o','*','s'};
markerStyles2 = {'s','o','*'};

Tcount = 16;
Tset = linspace(5,20,Tcount);

eps = 0.3;

%Hidden Markov Model (HMM)
eps_G_HMM = 0.07;
eps_B_HMM = 0.7;

%Markov model
eps_G_M = 0;
eps_B_M = 1;

NACK = 1;

rcount = 3;
rset = [0.1,0.2,0.3];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DELAY vs timeout T

meanDelay1 = zeros(Tcount,2*rcount);
Throughput1 = zeros(Tcount,2*rcount);


figure
clear str;      str = cell(1,2*rcount);
for i = 1:rcount
    r = rset(i);
    
    for t = 1:Tcount
        T = Tset(t);
        
        %HMM errors
        [P00,P01,P10,P11,P0x,Px0,P1x,Px1,P_kron,pi_kron,pi_I_kron] = TransitionProbabilityMatrices(eps_G_HMM,eps_B_HMM,eps,eps_G_HMM,eps_B_HMM,eps,rf,rr,NACK,K);                                                                                                   
        [PhiD1, phiD1, meanDelay1(t,2*i-1), varDelay1, Throughput1(t,2*i-1)] = NoCodingPhiD(z,k,K,T,P_kron,pi_I_kron,P00,P01,P10,P11,Px0,Px1,symbolic); 
        
        %Markov errors
        [P00,P01,P10,P11,P0x,Px0,P1x,Px1,P_kron,pi_kron,pi_I_kron] = TransitionProbabilityMatrices(eps_G_M,eps_B_M,eps,eps_G_M,eps_B_M,eps,rf,rr,NACK,K);                                                                                                  
        [PhiD1, phiD1, meanDelay1(t,2*i), varDelay1, Throughput1(t,2*i)] = NoCodingPhiD(z,k,K,T,P_kron,pi_I_kron,P00,P01,P10,P11,Px0,Px1,symbolic); 
        
    end  
end

%HMM errors
for i = 1:rcount
    str{i} = ['HMM, Uncoded ARQ, r=' num2str(rset(i))]; 
    plot(Tset,meanDelay1(:,2*i-1),'k','linewidth',2,...
        'color',colorStyles{1+rem(2*i-1,numel(colorStyles))},'linestyle','--',...
        'marker',markerStyles{1+rem(2*i-1,numel(markerStyles))},'markersize',8); hold on;
end
%Markov errors
for i = 1:rcount
    str{3+i} = ['Markov, Uncoded ARQ, r=' num2str(rset(i))]; 
    plot(Tset,meanDelay1(:,2*i),'k','linewidth',2,...
        'color',colorStyles{1+rem(2*i,numel(colorStyles))},'linestyle',':',...
        'marker',markerStyles2{1+rem(2*i,numel(markerStyles2))},'markersize',8); hold on; 
end
xlab = 'Timeout, $T$'; 
ylab = 'Average delay, $\bar{D}$';
box on;     set(gca,'FontSize',20) 
xlhand = get(gca,'xlabel'); xlabel(xlab,'Interpreter','latex'); set(xlhand,'fontsize',20) 
ylhand = get(gca,'ylabel'); ylabel(ylab,'Interpreter','latex'); set(ylhand,'fontsize',20)
legend(str,'FontSize',20,'location','NorthEast');
xaxis = 10;    yaxis = 8;
text(xaxis,yaxis,['\epsilon=' num2str(eps) ', \epsilon_G=' num2str(eps_Gf) ', \epsilon_B=' num2str(eps_Bf) ', k=' num2str(k)],'fontsize',20);
set(gca,'XLim',[5 20])
set(gca,'XTick',(5:5:20))
set(gca,'YLim',[7.5 12.5])
set(gca,'YTick',(7.5:0.5:12.5))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%THROUGHPUT vs timeout T

%HMM errors
figure
for i = 1:rcount
    plot(Tset,Throughput1(:,2*i-1),'k','linewidth',2,...
        'color',colorStyles{1+rem(2*i-1,numel(colorStyles))},'linestyle','--',...
        'marker',markerStyles{1+rem(2*i-1,numel(markerStyles))},'markersize',8); hold on; 
end
%Markov errors
for i = 1:rcount
    plot(Tset,Throughput1(:,2*i),'k','linewidth',2,...
        'color',colorStyles{1+rem(2*i,numel(colorStyles))},'linestyle',':',...
        'marker',markerStyles2{1+rem(2*i,numel(markerStyles2))},'markersize',8); hold on; 
end
xlab = 'Timeout, $T$'; 
ylab = 'Throughput, $\eta$';
box on;     set(gca,'FontSize',20) 
xlhand = get(gca,'xlabel'); xlabel(xlab,'Interpreter','latex'); set(xlhand,'fontsize',20) 
ylhand = get(gca,'ylabel'); ylabel(ylab,'Interpreter','latex'); set(ylhand,'fontsize',20)
legend(str,'FontSize',20,'location','NorthEast');
xaxis = 10;    yaxis = .52;
text(xaxis,yaxis,['\epsilon=' num2str(eps) ', \epsilon_G=' num2str(eps_Gf) ', \epsilon_B=' num2str(eps_Bf) ', k=' num2str(k)],'fontsize',20);
set(gca,'XLim',[5 20])
set(gca,'XTick',(5:5:20))
set(gca,'YLim',[0.5 0.71])
set(gca,'YTick',(0.5:0.02:0.71))
hold on; plot(Tset,(1-eps)*ones(1,Tcount),'k.')


%% 5-Coded ARQ scheme with HARQ combining
%This is left as future work.


%meanDelay4 = zeros(T_count,eps_count);
%Throughput4 = zeros(T_count,eps_count);
%[PhiD4, phiD4, meanDelay4(t,i), Throughput4(t,i)] = CodingHARQPhiD(...);








