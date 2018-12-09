%Another way to compute the average delay and the throughput for memoryless
%channels

clear; clc; 

lineStyles = {'--','-','--','-.',':'};
markerStyles = {'none'};
colorStyles = {'k'};
pos = [950,200,500,380];

eps_count = 100;
eps_set = linspace(0,0.5,eps_count) ; %the probability of block error (design requirement)

T_count = 3;
T_set = [15,10,8]; %timer=time to timeout (varies between k and 20)

k = 5;

delay1 = zeros(T_count,eps_count);
throughput1 = zeros(T_count,eps_count);
delay2 = zeros(T_count,eps_count);
throughput2 = zeros(T_count,eps_count);

for t = 1:T_count
    T = T_set(t);
    
    count = 0;
    for eps = eps_set
        count = count+1;
        
        d = T-k; %time to timer expiration
        
        
        %M=1 packet %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        delay1(t,count) = k-1+eps/(1-eps)*(k+(T-k)*eps)+1+eps*(1/(1-eps)+1);
        
        transmissiontime1 = (1/(1-eps) + eps^(d+1)*(1-eps^T))/(1-eps^(d+1+T));%forward link and reverse link                                                   
        %transmissiontime1 = 1/(1-eps)+eps^(d+1)*(eps^T/(1-eps^T)+1-eps^T);%1/((1-eps^(d+1))*(1-eps*(1-eps)));
        transmissiontime1 = 1/(1-eps)+eps^(d+1)*(1-eps^T)^(-1);     
        
        
        throughput1(t,count) = 1/transmissiontime1;
        
        
        %M=2 packets %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
        delayB1 = (k+(T-k)*eps)/(1-eps);
        delayB2 = (k+1+(T-k)*eps)/(1-eps^2);
        ackdelay = 1+eps*(1/(1-eps)+1);
        
        PA1 = 2*eps*(1-eps)^2; 
        PG3 = 2*eps^2*(1-eps);
        PA3A1 = PG3*(1+1/(1-eps^2))*eps*(1-eps);
       
        
        delay2(t,count) = 1+(k-1)...                                    %till A
                        + eps^2*delayB2...                              %self-loop B2
                        + PA1*(delayB1+ackdelay)...                     %A1 towards O
                        + PG3*( T+eps^2*(T+1)/(1-eps^2)...              %triangular loop (G3-B3-A3)
                        + eps*(1-eps)*(delayB1+ackdelay) ...            %G3, A3 and A1 towards 0
                        + (1-eps*(1-eps)-1/(1-eps^2))*ackdelay )...     %G3 and A3 towards O
                        + (1-PA1-PG3-PA3A1)*ackdelay;                   %to C2 and from C2 O
                                                         
        transmissiontime2 = 2+eps^2/(1-eps^2)...
                          + 2*eps*(1-eps)^2*(eps/(1-eps)+eps*(1-eps)*eps^(T-k)*(1/(1-eps^T)+1-eps^T))...
                          + eps*(1-eps)^2*eps^(T-k)*(1/(1-eps^T)+1-eps^T)...
                          + 2*eps^2*(1-eps)*(1+eps^2/(1-eps^2))...
                          *(eps*(1-eps)*(eps/(1-eps)+eps*(1-eps)*eps^d*(1/(1-eps^T)+1-eps^T))+eps*(1-eps)*eps^d*(1/(1-eps^T)+1-eps^T));

        throughput2(t,count) = 2/transmissiontime2;
    end
end


clear str;      str = cell(1,2*T_count);
   
%DELAY
figure
for t = 1:T_count
    T = T_set(t);
    
    str{2*t-1} = ['Uncoded ARQ, T=' num2str(T)]; 
    str{2*t} = ['Coded ARQ, T=' num2str(T)]; 

     
    plot(eps_set,delay1(t,:),'k','linewidth',2,...
        'color',colorStyles{1+rem(t,numel(colorStyles))},'linestyle',lineStyles{1+rem(t,numel(lineStyles))},...
        'marker',markerStyles{1+rem(t,numel(markerStyles))},'markersize',12); hold on;
    
    plot(eps_set,delay2(t,:),'k','linewidth',2,...
        'color','r','linestyle',lineStyles{1+rem(t,numel(lineStyles))},...
        'marker',markerStyles{1+rem(t,numel(markerStyles))},'markersize',12); hold on;
end
xlab = 'Erasure rate, $\epsilon$'; 
ylab = 'Average delay, $\bar{D}$';
box on;     set(gca,'FontSize',20) 
xlhand = get(gca,'xlabel'); xlabel(xlab,'Interpreter','latex'); set(xlhand,'fontsize',20) 
ylhand = get(gca,'ylabel'); ylabel(ylab,'Interpreter','latex'); set(ylhand,'fontsize',20)
legend(str,'FontSize',20,'location','NorthEast');
xaxis = 0.3;    yaxis = 2;
%text(xaxis,yaxis,['M=' num2str(M) ', N=' num2str(N) ', r=' num2str(r) ', k=' num2str(k)],'fontsize',20);


%THROUGHPUT
figure
for t = 1:T_count
    T = T_set(t);
     
    plot(eps_set,throughput1(t,:),'k','linewidth',2,...
        'color',colorStyles{1+rem(t,numel(colorStyles))},'linestyle',lineStyles{1+rem(t,numel(lineStyles))},...
        'marker',markerStyles{1+rem(t,numel(markerStyles))},'markersize',12); hold on;
    
    plot(eps_set,throughput2(t,:),'k','linewidth',2,...
        'color','r','linestyle',lineStyles{1+rem(t,numel(lineStyles))},...
        'marker',markerStyles{1+rem(t,numel(markerStyles))},'markersize',12); hold on;
end
xlab = 'Erasure rate, $\epsilon$'; 
ylab = 'Throughput, $\eta$';
box on;     set(gca,'FontSize',20) 
xlhand = get(gca,'xlabel'); xlabel(xlab,'Interpreter','latex'); set(xlhand,'fontsize',20) 
ylhand = get(gca,'ylabel'); ylabel(ylab,'Interpreter','latex'); set(ylhand,'fontsize',20)
legend(str,'FontSize',20,'location','NorthEast');
xaxis = 0.3;    yaxis = 2;
%text(xaxis,yaxis,['M=' num2str(M) ', N=' num2str(N) ', r=' num2str(r) ', k=' num2str(k)],'fontsize',20);

