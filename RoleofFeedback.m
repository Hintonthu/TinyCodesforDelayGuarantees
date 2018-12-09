%by Derya Malak
%Last revised: April 15, 2018 Sunday
%Role of Feedback in ARQ

clear; clc; %close all; 

lineStyles = {'--','-','--','-.',':'};
markerStyles = {'none'};
colorStyles = {'k'};
pos = [950,200,500,380];

%Transmission time (\tau): # of frames being transmitted per successful frame

%Total delay (D): time from when a frame is first transmitted to when its ACK received

%# of channel states of a multistate Markov process S_t
K = 2;      %There are Good G and Bad B states

rf = 0.1;    %1/r represents the average error burst    %(or r = 0.3)    
rr = 0.1;
r = rf;
if r == 0
    K = 1;
end


scale = 1;
k = scale*5;      %round trip time (RTT)


eps_count = 20;
eps_set = linspace(0.001,0.5,eps_count) ; %the probability of block error (design requirement)
z = 1;     %This is the z-transform parameter

T_count = 2;
T_set = scale*[15,8]; %timer=time to timeout (varies between k and 20)
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


%% 1-Uncoded ARQ scheme with no HARQ combining

meanDelay1 = zeros(T_count,eps_count);
Throughput1 = zeros(T_count,eps_count);

varDelay1 = zeros(T_count,eps_count);

%DELAY vs block-error rate \epsilon
figure
clear str1;      str1 = cell(1,T_count);
for t = 1:T_count
    T = T_set(t);
    str1{t} = ['Uncoded ARQ, T=' num2str(T)]; 
    
    for i = 1:eps_count
        epsf = eps_set(i);
        epsr = epsf;

        [P00,P01,P10,P11,P0x,Px0,P1x,Px1,P_kron,pi_kron,pi_I_kron] = TransitionProbabilityMatrices(eps_Gf,eps_Bf,epsf,eps_Gr,eps_Br,epsr,rf,rr,NACK);
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
legend(str1,'FontSize',20,'location','NorthEast');
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
legend(str1,'FontSize',20,'location','NorthEast');
xaxis = 0.3;    yaxis = 0.35;
text(xaxis,yaxis,['r=' num2str(r) ', k=' num2str(k)],'fontsize',20);


%PREDICTABLE DELAY vs block-error rate \epsilon
figure
for t = 1:T_count
    T = T_set(t);
    
    plot(eps_set,meanDelay1(t,:)+2*sqrt(varDelay1(t,:)),'k','linewidth',2,...
        'color',colorStyles{1+rem(t,numel(colorStyles))},'linestyle',lineStyles{1+rem(t,numel(lineStyles))},...
        'marker',markerStyles{1+rem(t,numel(markerStyles))},'markersize',12); hold on;
end
xlab = 'Erasure rate, $\epsilon$'; 
ylab = 'Predictable delay, $\hat{D}$';
box on;     set(gca,'FontSize',20) 
xlhand = get(gca,'xlabel'); xlabel(xlab,'Interpreter','latex'); set(xlhand,'fontsize',20) 
ylhand = get(gca,'ylabel'); ylabel(ylab,'Interpreter','latex'); set(ylhand,'fontsize',20)
legend(str1,'FontSize',20,'location','NorthEast');
xaxis = 0.3;    yaxis = 5;
text(xaxis,yaxis,['r=' num2str(r) ', k=' num2str(k)],'fontsize',20);




%% 3-Cumulative feedback-based ARQ (CF ARQ) scheme with no HARQ combining 

meanDelay3 = zeros(T_count,eps_count);
Throughput3 = zeros(T_count,eps_count);


%DELAY (average per packet) vs block-error rate \epsilon
figure
clear str2;      str2 = cell(1,T_count);
for t = 1:T_count
    T = T_set(t);
    str2{t} = ['CF ARQ, T=' num2str(T)]; 
    
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
          
        [P00,P01,P10,P11,P0x,Px0,P1x,Px1,P_kron,pi_kron,pi_I_kron] = TransitionProbabilityMatrices(eps_Gf,eps_Bf,epsf,eps_Gr,eps_Br,epsr,rf,rr,NACK);
        [PhiD3, phiD3, meanDelay3(t,i), Throughput3(t,i)] = CodingPhiD(z,k,K,T,P_kron,pi_I_kron,P00,P01,P10,P11,Px0,Px1,M,N,symbolic);

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
legend(str2,'FontSize',20,'location','NorthEast');
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
legend(str2,'FontSize',20,'location','NorthEast');
xaxis = 0.3;    yaxis = 2;
text(xaxis,yaxis,['M=' num2str(M) ', N=' num2str(N) ', r=' num2str(r) ', k=' num2str(k)],'fontsize',20);


%% ROLE OF FEEDBACK

eps = 0.1;
asymmetric = 1;  % Fixed erasure rate (eps_f=eps_r=eps)
if asymmetric ==1
    epsf = 0.5;
    epsr = 0.1;
else
    epsf = eps;
    epsr = eps;
end

eps_Gf = 0;  %no error in good state
eps_Bf = 1;  %always error in bad state

eps_Gr = 0;  %no error in good state
eps_Br = 1;  %always error in bad state

%Keep the forward (reverse) erasure rate same and increase the reverse (forward) erasure rate 
forward = 1;
r = 0.3;
if forward == 1
    rf = r;    %1/r represents the average error burst    %(or r = 0.3)    
else
    rr = r;
end

r_count = 99;
r_set = linspace(0.01,0.99,r_count);


% 1-Uncoded ARQ scheme with no HARQ combining

meanDelay1 = zeros(T_count,eps_count);
Throughput1 = zeros(T_count,eps_count);

varDelay1 = zeros(T_count,eps_count);

%DELAY vs block-error rate \epsilon
figure
clear str1;      str1 = cell(1,T_count);
for t = 1:T_count
    T = T_set(t);
    str1{t} = ['ARQ, T=' num2str(T)]; 
    
    for i = 1:r_count
        if forward == 1
            rr = r_set(i); 
            textstr1 = ['\epsilon=' num2str(eps) ', k=' num2str(k) ', r^{(f)}=' num2str(rf)];
            xlab = 'Average error burst of forward link, $1/r^{(r)}$'; 
        else
            rf = r_set(i); 
            textstr1 = ['\epsilon=' num2str(eps) ', k=' num2str(k) ', r^{(r)}=' num2str(rr)];
            xlab = 'Average error burst of forward link, $1/r^{(f)}$'; 
        end

        [P00,P01,P10,P11,P0x,Px0,P1x,Px1,P_kron,pi_kron,pi_I_kron] = TransitionProbabilityMatrices(eps_Gf,eps_Bf,epsf,eps_Gr,eps_Br,epsr,rf,rr,NACK);
        [PhiD1, phiD1, meanDelay1(t,i), varDelay1(t,i), Throughput1(t,i)] = NoCodingPhiD(z,k,K,T,P_kron,pi_I_kron,P00,P01,P10,P11,Px0,Px1,symbolic);  
    end
    plot(1./r_set,meanDelay1(t,:),'k','linewidth',2,...
        'linestyle',':',...
        'marker',markerStyles{1+rem(t,numel(markerStyles))},'markersize',12); hold on;
end 
ylab = 'Average delay, $\bar{D}$';
box on;     set(gca,'FontSize',20) 
xlhand = get(gca,'xlabel'); xlabel(xlab,'Interpreter','latex'); set(xlhand,'fontsize',20) 
ylhand = get(gca,'ylabel'); ylabel(ylab,'Interpreter','latex'); set(ylhand,'fontsize',20)
legend(str1,'FontSize',20,'location','NorthEast');
xaxis = 10;    yaxis = 10;
text(xaxis,yaxis,textstr1,'fontsize',20);



%THROUGHPUT vs block-error rate \epsilon
figure
for t = 1:T_count
    T = T_set(t);
     
    plot(1./r_set,Throughput1(t,:),'k','linewidth',2,...
        'linestyle',':',...
        'marker',markerStyles{1+rem(t,numel(markerStyles))},'markersize',12); hold on;
end
ylab = 'Throughput, $\eta$';
box on;     set(gca,'FontSize',20) 
xlhand = get(gca,'xlabel'); xlabel(xlab,'Interpreter','latex'); set(xlhand,'fontsize',20) 
ylhand = get(gca,'ylabel'); ylabel(ylab,'Interpreter','latex'); set(ylhand,'fontsize',20)
legend(str1,'FontSize',20,'location','NorthEast');
xaxis = 10;    yaxis = 0.7;
text(xaxis,yaxis,textstr1,'fontsize',20);


%PREDICTABLE DELAY vs block-error rate \epsilon
figure
for t = 1:T_count
    T = T_set(t);
    
    plot(1./r_set,meanDelay1(t,:)+2*sqrt(varDelay1(t,:)),'k','linewidth',2,...
        'linestyle',':',...
        'marker',markerStyles{1+rem(t,numel(markerStyles))},'markersize',12); hold on;
end 
ylab = 'Predictable delay, $\hat{D}$';
box on;     set(gca,'FontSize',20) 
xlhand = get(gca,'xlabel'); xlabel(xlab,'Interpreter','latex'); set(xlhand,'fontsize',20) 
ylhand = get(gca,'ylabel'); ylabel(ylab,'Interpreter','latex'); set(ylhand,'fontsize',20)
legend(str1,'FontSize',20,'location','NorthEast');
xaxis = 10;    yaxis = 45;
text(xaxis,yaxis,textstr1,'fontsize',20);


% 3-Cumulative feedback-based ARQ (CF ARQ) scheme with no HARQ combining 

meanDelay3 = zeros(T_count,eps_count);
Throughput3 = zeros(T_count,eps_count);


%These values are modified for the CF model
epsf = sqrt(epsf^4+2*epsf^3*(1-epsf)); 
epsr = sqrt(epsr^4+2*epsr^3*(1-epsr)); 


%DELAY (average per packet) vs block-error rate \epsilon
figure
clear str2;      str2 = cell(1,T_count);
for t = 1:T_count
    T = T_set(t);
    str2{t} = ['CF ARQ, T=' num2str(T)]; 
    
    for i = 1:r_count       
        if forward == 1
            rr = r_set(i); 
            textstr2 = ['M=' num2str(M) ', N=' num2str(N) ', \epsilon=' num2str(eps) ', k=' num2str(k) ', r^{(f)}=' num2str(rf)];     
            xlab = 'Average error burst of feedback, $1/r^{(r)}$'; 
        else
            rf = r_set(i); 
            textstr2 = ['M=' num2str(M) ', N=' num2str(N) ', \epsilon=' num2str(eps) ', k=' num2str(k) ', r^{(r)}=' num2str(rr)];   
            xlab = 'Average error burst of forward link, $1/r^{(f)}$'; 
        end
                        
        [P00,P01,P10,P11,P0x,Px0,P1x,Px1,P_kron,pi_kron,pi_I_kron] = TransitionProbabilityMatrices(eps_Gf,eps_Bf,epsf,eps_Gr,eps_Br,epsr,rf,rr,NACK);
        [PhiD3, phiD3, meanDelay3(t,i), Throughput3(t,i)] = CodingPhiD(z,k,K,T,P_kron,pi_I_kron,P00,P01,P10,P11,Px0,Px1,M,N,symbolic);

    end
    plot(1./r_set,meanDelay3(t,:),'m','linewidth',2,...
        'linestyle','-.',...
        'marker',markerStyles{1+rem(t,numel(markerStyles))},'markersize',12); hold on;
end
ylab = 'Average delay, $\bar{D}$';
box on;     set(gca,'FontSize',20) 
xlhand = get(gca,'xlabel'); xlabel(xlab,'Interpreter','latex'); set(xlhand,'fontsize',20) 
ylhand = get(gca,'ylabel'); ylabel(ylab,'Interpreter','latex'); set(ylhand,'fontsize',20)
legend(str2,'FontSize',20,'location','NorthEast');
xaxis = 10;    yaxis = 10;
text(xaxis,yaxis,textstr2,'fontsize',20);

 
%THROUGHPUT vs block-error rate \epsilon
figure
for t = 1:T_count
    T = T_set(t);
      
     plot(1./r_set,Throughput3(t,:),'m','linewidth',2,...
        'linestyle','-.',...
        'marker',markerStyles{1+rem(t,numel(markerStyles))},'markersize',12); hold on;
end 
ylab = 'Throughput, $\eta$';
box on;     set(gca,'FontSize',20) 
xlhand = get(gca,'xlabel'); xlabel(xlab,'Interpreter','latex'); set(xlhand,'fontsize',20) 
ylhand = get(gca,'ylabel'); ylabel(ylab,'Interpreter','latex'); set(ylhand,'fontsize',20)
legend(str2,'FontSize',20,'location','NorthEast');
xaxis = 10;    yaxis = 0.7;
text(xaxis,yaxis,textstr2,'fontsize',20);



%uncoded: as r decreases delay increases (mostly), CF opposite
%uncoded: as r decreases throughput decreases, CF opposite (but more stable)

%% How to compensate the forward errors with feedback? 
%How much forward erasures can be tolerated with CF?

eps = 0.1;

eps_Gf = 0;  %no error in good state
eps_Bf = 1;  %always error in bad state

eps_Gr = 0;  %no error in good state
eps_Br = 1;  %always error in bad state

%Keep the reverse erasure rate same and increase the reverse (forward) erasure rate 
forward = 0;
rr = 0.1;
rf = 0.3;
epsr = eps;


eps_count = 20;
eps_set = linspace(0.001,0.5,eps_count) ; %the probability of block error (design requirement)



% 1-Uncoded ARQ scheme with no HARQ combining

meanDelay1 = zeros(T_count,eps_count);
Throughput1 = zeros(T_count,eps_count);

varDelay1 = zeros(T_count,eps_count);

%DELAY vs block-error rate \epsilon
figure
clear str1;      str1 = cell(1,T_count);
for t = 1:T_count
    T = T_set(t);
    str1{t} = ['ARQ, T=' num2str(T)]; 
    
    for i = 1:eps_count
        if forward == 1
            rr = r_set(i); 
            textstr1 = ['\epsilon=' num2str(eps) ', k=' num2str(k) ', r^{(f)}=' num2str(rf)];
            xlab = 'Average error burst of forward link, $1/r^{(r)}$'; 
        else
            epsf = eps_set(i); 
            textstr1 = ['\epsilon^{(r)}=' num2str(epsr) ', k=' num2str(k) ', r^{(r)}=' num2str(rr) ', r^{(f)}=' num2str(rf)];
            xlab = 'Erasure rate of forward link, $\epsilon^{(f)}$';
        end

        [P00,P01,P10,P11,P0x,Px0,P1x,Px1,P_kron,pi_kron,pi_I_kron] = TransitionProbabilityMatrices(eps_Gf,eps_Bf,epsf,eps_Gr,eps_Br,epsr,rf,rr,NACK);
        [PhiD1, phiD1, meanDelay1(t,i), varDelay1(t,i), Throughput1(t,i)] = NoCodingPhiD(z,k,K,T,P_kron,pi_I_kron,P00,P01,P10,P11,Px0,Px1,symbolic);  
    end
    plot(eps_set,meanDelay1(t,:),'k','linewidth',2,...
        'linestyle',':',...
        'marker',markerStyles{1+rem(t,numel(markerStyles))},'markersize',12); hold on;
end 
ylab = 'Average delay, $\bar{D}$';
box on;     set(gca,'FontSize',20) 
xlhand = get(gca,'xlabel'); xlabel(xlab,'Interpreter','latex'); set(xlhand,'fontsize',20) 
ylhand = get(gca,'ylabel'); ylabel(ylab,'Interpreter','latex'); set(ylhand,'fontsize',20)
legend(str1,'FontSize',20,'location','NorthEast');
xaxis = .1;    yaxis = 8;
text(xaxis,yaxis,textstr1,'fontsize',20);



%THROUGHPUT vs block-error rate \epsilon
figure
for t = 1:T_count
    T = T_set(t);
     
    plot(eps_set,Throughput1(t,:),'k','linewidth',2,...
        'linestyle',':',...
        'marker',markerStyles{1+rem(t,numel(markerStyles))},'markersize',12); hold on;
end
ylab = 'Throughput, $\eta$';
box on;     set(gca,'FontSize',20) 
xlhand = get(gca,'xlabel'); xlabel(xlab,'Interpreter','latex'); set(xlhand,'fontsize',20) 
ylhand = get(gca,'ylabel'); ylabel(ylab,'Interpreter','latex'); set(ylhand,'fontsize',20)
legend(str1,'FontSize',20,'location','NorthEast');
xaxis = 0.1;    yaxis = 0.7;
text(xaxis,yaxis,textstr1,'fontsize',20);


%PREDICTABLE DELAY vs block-error rate \epsilon
figure
for t = 1:T_count
    T = T_set(t);
    
    plot(eps_set,meanDelay1(t,:)+2*sqrt(varDelay1(t,:)),'k','linewidth',2,...
        'linestyle',':',...
        'marker',markerStyles{1+rem(t,numel(markerStyles))},'markersize',12); hold on;
end 
ylab = 'Predictable delay, $\hat{D}$';
box on;     set(gca,'FontSize',20) 
xlhand = get(gca,'xlabel'); xlabel(xlab,'Interpreter','latex'); set(xlhand,'fontsize',20) 
ylhand = get(gca,'ylabel'); ylabel(ylab,'Interpreter','latex'); set(ylhand,'fontsize',20)
legend(str1,'FontSize',20,'location','NorthEast');
xaxis = .1;    yaxis = 5;
text(xaxis,yaxis,textstr1,'fontsize',20);


% 3-Cumulative feedback-based ARQ (CF ARQ) scheme with no HARQ combining 

meanDelay3 = zeros(T_count,eps_count);
Throughput3 = zeros(T_count,eps_count);


%These values are modified for the CF model
epsr = sqrt(epsr^4+2*epsr^3*(1-epsr)); 


%DELAY (average per packet) vs block-error rate \epsilon
figure
clear str2;      str2 = cell(1,T_count);
for t = 1:T_count
    T = T_set(t);
    str2{t} = ['CF ARQ, T=' num2str(T)]; 
    
    for i = 1:eps_count       
        if forward == 1
            rr = r_set(i); 
            textstr2 = ['M=' num2str(M) ', N=' num2str(N) ', \epsilon=' num2str(eps) ', k=' num2str(k) ', r^{(f)}=' num2str(rf)];     
            xlab = 'Average error burst of feedback, $1/r^{(r)}$'; 
        else
            epsf = eps_set(i); 
            epsf = sqrt(epsf^4+2*epsf^3*(1-epsf)); 
            
            textstr2 = ['M=' num2str(M) ', N=' num2str(N) ', \epsilon^{(r)}=' num2str(epsr) ', k=' num2str(k) ', r^{(r)}=' num2str(rr) ', r^{(f)}=' num2str(rf)];   
            xlab = 'Erasure rate of forward link, $\epsilon^{(f)}$';
        end
                        
        [P00,P01,P10,P11,P0x,Px0,P1x,Px1,P_kron,pi_kron,pi_I_kron] = TransitionProbabilityMatrices(eps_Gf,eps_Bf,epsf,eps_Gr,eps_Br,epsr,rf,rr,NACK);
        [PhiD3, phiD3, meanDelay3(t,i), Throughput3(t,i)] = CodingPhiD(z,k,K,T,P_kron,pi_I_kron,P00,P01,P10,P11,Px0,Px1,M,N,symbolic);

    end
    plot(eps_set,meanDelay3(t,:),'m','linewidth',2,...
        'linestyle','-.',...
        'marker',markerStyles{1+rem(t,numel(markerStyles))},'markersize',12); hold on;
end
ylab = 'Average delay, $\bar{D}$';
box on;     set(gca,'FontSize',20) 
xlhand = get(gca,'xlabel'); xlabel(xlab,'Interpreter','latex'); set(xlhand,'fontsize',20) 
ylhand = get(gca,'ylabel'); ylabel(ylab,'Interpreter','latex'); set(ylhand,'fontsize',20)
legend(str2,'FontSize',20,'location','NorthEast');
xaxis = .1;    yaxis = 8;
text(xaxis,yaxis,textstr2,'fontsize',20);

 
%THROUGHPUT vs block-error rate \epsilon
figure
for t = 1:T_count
    T = T_set(t);
      
     plot(eps_set,Throughput3(t,:),'m','linewidth',2,...
        'linestyle','-.',...
        'marker',markerStyles{1+rem(t,numel(markerStyles))},'markersize',12); hold on;
end 
ylab = 'Throughput, $\eta$';
box on;     set(gca,'FontSize',20) 
xlhand = get(gca,'xlabel'); xlabel(xlab,'Interpreter','latex'); set(xlhand,'fontsize',20) 
ylhand = get(gca,'ylabel'); ylabel(ylab,'Interpreter','latex'); set(ylhand,'fontsize',20)
legend(str2,'FontSize',20,'location','NorthEast');
xaxis = .1;    yaxis = 0.5;
text(xaxis,yaxis,textstr2,'fontsize',20);

