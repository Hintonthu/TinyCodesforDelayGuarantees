%by Derya Malak
%Last revised: Dec 8, 2018 
%Simulation of the ARQ models

clear; clc; %close all


lineStyles = {'none'};
markerStyles = {'o'};
colorStyles = {'k'};
pos = [950,200,500,380];

cnst = 1;
k = cnst*5; %RTT
%T = cnst*8; %timeout value 1
T = cnst*15; %timeout value 2

d = T-k; %time difference between timeout T and RTT k

%Forward and reverse channel parameters (symmetric GE channels)
%r = 0.3;
r = 0.1;

eps_count = 20;
eps_min = 0.001;
eps_max = 0.5;
eps_set = linspace(eps_min,eps_max,eps_count) ; %the probability of block error (design requirement)


eps_Gf = 0;  %no error in good state
eps_Bf = 1;  %always error in bad state

eps_Gr = 0;  %no error in good state
eps_Br = 1;  %always error in bad state


%Coded ARQ Model
M = 2;      %total number of transmitted packets
N = 2;      %required number of successfully received packets   

NACK = 1; %1 with NACK and 0 without NACK

symbolic = 0; %1 for memoryless channel results

%NOTE THAT TO BE ABLE GET GOOD AVERAGED RESULTS, THE FOLLOWING NUMBERS
%SHOULD BE LARGE ENOUGH
packetno = 10^4;    %10^4;   %10^4; %ideally all of them are successfully transmitted when eps=0
slotno = 10^6;      %10^6;   %number of consecutive observations from the channel 
total_r_no = 10;     %100;    %total number of channel realizations

% if r == 0.1
%     slotno = 10*slotno;
% end
if T == 8
    total_r_no = total_r_no*10;
end


%NOTATION
%Transmission time tau: number of packets transmitted per successful packet 
TransmissionTime = zeros(eps_count,packetno);

%Throughput eta=1/E[tau]: reciprocal of the average transmission time
%Delay D: time from when a packet is first transmitted to when its ACK is successfully received at the sender
Delay = zeros(eps_count,packetno);



%STATE I: transmission of a new packet by the sender

%THROUGHPUT

%TRANSITION TO STATE A: after sending a new packet, the transmitter receives a fb k-1 time slots later. 

%TRANSITION TO STATE B: retransmission of an erroneous packet

%TRANSITION TO STATE C: retransmission of a packet that was correctly received, but with its timer expiring

%TRANSITION TO STATE O: if fb is error-free ACK  (received before timer expiration)

%DELAY

%TRANSITION TO STATE G: If the feedback is an erroneous NACK, the transmitter waits for timeout. 
%The packet will be retransmitted after the timer expires, which involves a delay (transition from G to B).


%% RUN THE FOLLOWING SECTIONS INDIVIDUALLY TO CHARACTERIZE THE PERFORMANCE OF THE DELAY AND THE THROUGHPUT
% GILBERT ELLIOTT CHANNEL (I.E. WHEN K=2)

% THE FOLLOWING FUNCTIONS RETURN THE AVERAGE DELAY, AND THROUGHPUT FOR THE RESPECTIVE SCHEMES

% 1- FOR UNCODED ARQ
% 2- FOR UNCODED HARQ
% 3- FOR CUMULATIVE FEEDBACK ARQ
% 4- FOR CODED ARQ
    
%% 1-Uncoded ARQ scheme with no HARQ combining

TransmissionTime = zeros(eps_count,packetno);

Delay = zeros(eps_count,packetno);

qset = eps_set*r./(1-eps_set);  %since eps_B q/(q+r) = eps, and eps_B = 1


count = 0;
for q = qset
    count = count + 1
    
    for r_no = 1:total_r_no
        %r_no
        
        %Evolution of GE channel for total_no time slots (Single Realization)
        [ChannelStatef,ChannelStater,Pkron] = GEchannelUncodedARQ(q,r,slotno);
        
        total_no = 10^5;         

        for p = 1:packetno %Each packet can be considered independently of the others            
            %p
            
            t_tau = 1; %time slot for throughput flow graph (need to transmit at least once)
            
            t_D = p; %time slot for delay flow graph=(initial) transmission time slot for the p th packet

            total_no = total_no + p;

            %STATE I
            ACKreceived = 0;  

            %STATE A
            t_D = t_D+k-1;            
                       
            %TRANSITION FROM STATE A TO STATE B (SELF LOOP AROUND A)
            while ChannelStatef(t_D) == 0 && max(t_D,t_tau) < 1/3*total_no %Retransmit packet while forward channel is BAD
                
                t_D = t_D+1;
                if ChannelStater(t_D) == 1 %NACK received (immediate retransmission)

                    t_D = t_D+k-1;

                %STATE G for the DELAY FLOW GRAPH    
                else %if max(ChannelStater(t_D+1:t_D+d)) == 0  %NACK is LOST (wait till timeout (d slots) and retransmit packet) 

                    t_D = t_D+d+k-1;
                end    
                t_tau = t_tau + 1; 
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           
            t_D_for_tau = t_D;
            %DELAY GRAPH
            %TRANSITION FROM STATE A TO STATE 0 
            if ChannelStater(t_D) == 1
                t_D = t_D+1;
            else    
                %TRANSITION FROM STATE A TO STATE C 
                while ChannelStater(t_D) == 0 && max(t_D,t_tau) < 2/3*total_no  %ACK is LOST                  
                    t_D = t_D+1;                   
                end
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %THROUGHPUT GRAPH
            while ACKreceived == 0 && max(t_D_for_tau,t_tau) < total_no  %packet will be retransmitted till the ACK is successfully received.
                
                %TRANSITION FROM STATE A TO STATE C THEN TO STATE 0
                if max(ChannelStater(t_tau:t_tau+d)) == 0 %ACK is LOST and the timer expires   
                    %t_tau = t_tau + 1; %not sure about this?
                    while max(ChannelStater(t_tau:t_tau+T-1)) == 0 && max(t_D_for_tau,t_tau) < total_no 
                        t_tau = t_tau + 1;  %wait till timeout (d=T-k slots) and retransmit packet
                    end
                    ACKreceived = 1;
                
                
                %TRANSITION FROM STATE A TO STATE 0 
                elseif max(ChannelStater(t_tau:t_tau+d)) == 1 && max(t_D_for_tau,t_tau) < total_no %ACK arrives before timer expires

                    ACKreceived = 1;
             
                end
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            DelayperRealization = t_D-p;
            TransmissionTimeperRealization = t_tau;
            
            Delay(count,p) = Delay(count,p)+DelayperRealization;
            TransmissionTime(count,p) = TransmissionTime(count,p)+TransmissionTimeperRealization;
        end
    end
            
    Delay(count,:) = Delay(count,:)/total_r_no;
    TransmissionTime(count,:) = TransmissionTime(count,:)/total_r_no;
end

MeanDelay = zeros(1,eps_count);
ThroughputInv = zeros(1,eps_count);
for i = 1:eps_count
    for j = 1:packetno                
        MeanDelay(i) = MeanDelay(i) + Delay(i,j);
        ThroughputInv(i) = ThroughputInv(i) + TransmissionTime(i,j);
    end
    MeanDelay(i) = MeanDelay(i)/packetno;
    ThroughputInv(i) = ThroughputInv(i)/packetno;
end
Throughput = 1./ThroughputInv;
    
t = 1;
T_count = 1;
clear str;  str = cell(1,T_count);
str{t} = ['ARQ, T=' num2str(T)]; 

figure    
plot(eps_set,MeanDelay,'k','linewidth',2,...
        'color',colorStyles{1+rem(t,numel(colorStyles))},'linestyle',lineStyles{1+rem(t,numel(lineStyles))},...
        'marker',markerStyles{1+rem(t,numel(markerStyles))},'markersize',10); hold on;  
xlab = 'Erasure rate, $\epsilon$'; 
ylab = 'Average delay, $\bar{D}$';
box on;     set(gca,'FontSize',20) 
xlhand = get(gca,'xlabel'); xlabel(xlab,'Interpreter','latex'); set(xlhand,'fontsize',20) 
ylhand = get(gca,'ylabel'); ylabel(ylab,'Interpreter','latex'); set(ylhand,'fontsize',20)
legend(str,'FontSize',20,'location','NorthEast');
xaxis = 0.3;    yaxis = 5;
text(xaxis,yaxis,['M=' num2str(M) ', N=' num2str(N) ', r=' num2str(r) ', k=' num2str(k)],'fontsize',20);
   
figure
plot(eps_set,Throughput,'k','linewidth',2,...
        'color',colorStyles{1+rem(t,numel(colorStyles))},'linestyle',lineStyles{1+rem(t,numel(lineStyles))},...
        'marker',markerStyles{1+rem(t,numel(markerStyles))},'markersize',10); hold on;    
xlab = 'Erasure rate, $\epsilon$'; 
ylab = 'Throughput, $\eta$';
box on;     set(gca,'FontSize',20) 
xlhand = get(gca,'xlabel'); xlabel(xlab,'Interpreter','latex'); set(xlhand,'fontsize',20) 
ylhand = get(gca,'ylabel'); ylabel(ylab,'Interpreter','latex'); set(ylhand,'fontsize',20)
legend(str,'FontSize',20,'location','NorthEast');
xaxis = 0.3;    yaxis = 0.35;
text(xaxis,yaxis,['M=' num2str(M) ', N=' num2str(N) ', r=' num2str(r) ', k=' num2str(k)],'fontsize',20);
   

%% 2-Uncoded ARQ scheme and HARQ with soft combining

TransmissionTime = zeros(eps_count,packetno);

Delay = zeros(eps_count,packetno);

a = @(eps) -log(1-eps); %This is to make sure that eps_B(1) = eps_B
eps_HARQ = @(m,eps) 1-exp(-a(eps)./m); %epsB as function of the retransmission attempt m

%{
count = 0;
for eps = eps_set
    count = count + 1
    
    for r_no = 1:total_r_no
        %r_no

        q_HARQ = @(m) r*eps_HARQ(m,eps)./(1-eps_HARQ(m,eps));  %since eps_B q/(q+r) = eps
        q = eps*r/(1-eps);        
    
        %Evolution of GE channel for total_no time slots (Single Realization)
        [ChannelStatef,ChannelStater,Pkron] = GEchannelUncodedHARQ(q,q_HARQ,r,slotno);
        
        total_no = 10^5;         

        for p = 1:packetno %Each packet can be considered independently of the others            
            %p
            
            t_tau = 1; %time slot for throughput flow graph (need to transmit at least once)
            
            t_D = p; %time slot for delay flow graph=(initial) transmission time slot for the p th packet

            total_no = total_no + p;

            %STATE I
            ACKreceived = 0;  

            %STATE A
            t_D = t_D+k-1;            
                       
            %TRANSITION FROM STATE A TO STATE B (SELF LOOP AROUND A)
            while ChannelStatef(t_D) == 0 && max(t_D,t_tau) < 1/3*total_no %Retransmit packet while forward channel is BAD
                
                t_D = t_D+1;
                if ChannelStater(t_D) == 1 %NACK received (immediate retransmission)

                    t_D = t_D+k-1;

                %STATE G for the DELAY FLOW GRAPH    
                else %if max(ChannelStater(t_D+1:t_D+d)) == 0  %NACK is LOST (wait till timeout (d slots) and retransmit packet) 

                    t_D = t_D+d+k-1;
                end    
                t_tau = t_tau + 1; 
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            t_D_for_tau = t_D;
            %DELAY GRAPH
            %TRANSITION FROM STATE A TO STATE 0 

            t_D = t_D+1;
            %TRANSITION FROM STATE A TO STATE C 
            while ChannelStater(t_D) == 0 && max(t_D,t_tau) < 2/3*total_no  %ACK is LOST                  
                t_D = t_D+1;                   
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %THROUGHPUT GRAPH
            while ACKreceived == 0 && max(t_D_for_tau,t_tau) < total_no  %packet will be retransmitted till the ACK is successfully received.
                
                %TRANSITION FROM STATE A TO STATE C THEN TO STATE 0
                if max(ChannelStater(t_tau:t_tau+d)) == 0 %ACK is LOST and the timer expires   
                    t_tau = t_tau + 1; %not sure about this?
                    while max(ChannelStater(t_tau:t_tau+T-1)) == 0 && max(t_D_for_tau,t_tau) < total_no 
                        t_tau = t_tau + 1;  %wait till timeout (d=T-k slots) and retransmit packet
                    end
                    ACKreceived = 1;
                
                
                %TRANSITION FROM STATE A TO STATE 0 
                elseif max(ChannelStater(t_tau:t_tau+d)) == 1 && max(t_D_for_tau,t_tau) < total_no %ACK arrives before timer expires

                    ACKreceived = 1;
             
                end
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            DelayperRealization = t_D-p;
            TransmissionTimeperRealization = t_tau;
            
            Delay(count,p) = Delay(count,p)+DelayperRealization;
            TransmissionTime(count,p) = TransmissionTime(count,p)+TransmissionTimeperRealization;
        end
    end
            
    Delay(count,:) = Delay(count,:)/total_r_no;
    TransmissionTime(count,:) = TransmissionTime(count,:)/total_r_no;
    
end
%}

count = 0;
for eps = eps_set
    count = count + 1
    
    for r_no = 1:total_r_no
        %r_no

        q_HARQ = @(m) r*eps_HARQ(m,eps)./(1-eps_HARQ(m,eps));  %since eps_B q/(q+r) = eps, and eps_B = 1
        q = eps*r/(1-eps);        
        
        %Evolution of GE channel for total_no time slots (Single Realization)
        [ChannelStatef,ChannelStater,Pkron] = GEchannelUncodedHARQ(q,q_HARQ,r,slotno);
        
        total_no = 10^5;         

        for p = 1:packetno %Each packet can be considered independently of the others            
            %p                    

            t_tau = 1; %time slot for throughput flow graph (need to transmit at least once)
            
            t_D = p; %time slot for delay flow graph=(initial) transmission time slot for the p th packet

            total_no = total_no + p;

            %STATE I
            ACKreceived = 0;  

            %STATE A
            t_D = t_D+k-1;            
                       
            %TRANSITION FROM STATE A TO STATE B (SELF LOOP AROUND A)
            while ChannelStatef(t_D) == 0 && max(t_D,t_tau) < 1/3*total_no %Retransmit packet while forward channel is BAD
                               
                t_D = t_D+1;
                if ChannelStater(t_D) == 1 %NACK received (immediate retransmission)

                    t_D = t_D+k-1;

                %STATE G for the DELAY FLOW GRAPH    
                else %if max(ChannelStater(t_D+1:t_D+d)) == 0  %NACK is LOST (wait till timeout (d slots) and retransmit packet) 

                    t_D = t_D+d+k-1;
                end    
                t_tau = t_tau + 1; 
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            t_D_for_tau = t_D;
            %DELAY GRAPH
            %TRANSITION FROM STATE A TO STATE 0 

            if ChannelStater(t_D) == 1
                t_D = t_D+1;
            %TRANSITION FROM STATE A TO STATE C 
            else
                while ChannelStater(t_D) == 0 && max(t_D,t_tau) < 2/3*total_no  %ACK is LOST                  
                    t_D = t_D+1;                   
                end
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %THROUGHPUT GRAPH
            while ACKreceived == 0 && max(t_D_for_tau,t_tau) < total_no  %packet will be retransmitted till the ACK is successfully received.
                
                %TRANSITION FROM STATE A TO STATE C THEN TO STATE 0
                if max(ChannelStater(t_tau:t_tau+d)) == 0 %ACK is LOST and the timer expires   
                    t_tau = t_tau + 1; %not sure about this?
                    while max(ChannelStater(t_tau:t_tau+T-1)) == 0 && max(t_D_for_tau,t_tau) < total_no 
                        t_tau = t_tau + 1;  %wait till timeout (d=T-k slots) and retransmit packet
                    end
                    ACKreceived = 1;
                
                
                %TRANSITION FROM STATE A TO STATE 0 
                elseif max(ChannelStater(t_tau:t_tau+d)) == 1 && max(t_D_for_tau,t_tau) < total_no %ACK arrives before timer expires

                    ACKreceived = 1;
             
                end
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            DelayperRealization = t_D-p;
            TransmissionTimeperRealization = t_tau;
            
            Delay(count,p) = Delay(count,p)+DelayperRealization;
            TransmissionTime(count,p) = TransmissionTime(count,p)+TransmissionTimeperRealization;
        end
    end
            
    Delay(count,:) = Delay(count,:)/total_r_no;
    TransmissionTime(count,:) = TransmissionTime(count,:)/total_r_no;
    
end

MeanDelay = zeros(1,eps_count);
ThroughputInv = zeros(1,eps_count);
for i = 1:eps_count
    for j = 1:packetno                
        MeanDelay(i) = MeanDelay(i) + Delay(i,j);
        ThroughputInv(i) = ThroughputInv(i) + TransmissionTime(i,j);
    end
    MeanDelay(i) = MeanDelay(i)/packetno;
    ThroughputInv(i) = ThroughputInv(i)/packetno;
end
Throughput = 1./ThroughputInv;
    
t = 1;
T_count = 1;
clear str;  str = cell(1,T_count);
str{t} = ['HARQ, T=' num2str(T)]; 

figure    
plot(eps_set,MeanDelay,'k','linewidth',2,...
        'color',colorStyles{1+rem(t,numel(colorStyles))},'linestyle',lineStyles{1+rem(t,numel(lineStyles))},...
        'marker',markerStyles{1+rem(t,numel(markerStyles))},'markersize',10); hold on;  
xlab = 'Erasure rate, $\epsilon$'; 
ylab = 'Average delay, $\bar{D}$';
box on;     set(gca,'FontSize',20) 
xlhand = get(gca,'xlabel'); xlabel(xlab,'Interpreter','latex'); set(xlhand,'fontsize',20) 
ylhand = get(gca,'ylabel'); ylabel(ylab,'Interpreter','latex'); set(ylhand,'fontsize',20)
legend(str,'FontSize',20,'location','NorthEast');
xaxis = 0.3;    yaxis = 5;
text(xaxis,yaxis,['M=' num2str(M) ', N=' num2str(N) ', r=' num2str(r) ', k=' num2str(k)],'fontsize',20);
   
figure
plot(eps_set,Throughput,'k','linewidth',2,...
        'color',colorStyles{1+rem(t,numel(colorStyles))},'linestyle',lineStyles{1+rem(t,numel(lineStyles))},...
        'marker',markerStyles{1+rem(t,numel(markerStyles))},'markersize',10); hold on;    
xlab = 'Erasure rate, $\epsilon$'; 
ylab = 'Throughput, $\eta$';
box on;     set(gca,'FontSize',20) 
xlhand = get(gca,'xlabel'); xlabel(xlab,'Interpreter','latex'); set(xlhand,'fontsize',20) 
ylhand = get(gca,'ylabel'); ylabel(ylab,'Interpreter','latex'); set(ylhand,'fontsize',20)
legend(str,'FontSize',20,'location','NorthEast');
xaxis = 0.3;    yaxis = 0.35;
text(xaxis,yaxis,['M=' num2str(M) ', N=' num2str(N) ', r=' num2str(r) ', k=' num2str(k)],'fontsize',20);


%% 3-Cumulative feedback-based ARQ (CF ARQ) scheme with no HARQ combining

%M=2 packets and N=2

%Combined channel probabilities
%{
%BOTH PACKETS NEED TO BE RETRANSMITTED
PCF10(1) = P10 * P10 + P10 * P01 + P01 * P10;  %FROM STATE A2 TO STATE B2
%probability: 3*eps^2*(1-eps)^2
PCF11(1) = P11 * P11 + P11 * P10 + P10 * P11;  %FROM STATE A2 TO STATE G2
%probability: eps^4+2*eps^3*(1-eps) 
eps_set_CF_f = eps_set.^4+2*eps_set.^3.*(1-eps_set) + 3*eps_set.^2.*(1-eps_set).^2;

P1xC = P10C+P11C;

%BOTH PACKETS ARE SUCCESSFULLY TRANSMITTED
PCF00(1) = P00 * P00;                          %FROM STATE A2 TO STATE 0
%probability: (1-eps)^4
PCF01(1) = P01 * P01 + P01 * P00 + P00 * P01;  %FROM STATE A2 TO STATE C2
%probability: (1-eps)^2*eps^2+2*(1-eps)^3*eps
eps_set_CF_r = eps_set.^4+2*eps_set.^3.*(1-eps_set) + (1-eps_set).^2.*eps_set.^2 + 2*(1-eps_set).^3.*eps_set;

P0xC = P00C+P01C;

%ONLY EITHER 1 OF THE PACKETS IS CORRECTLY ACKNOWLEDGED
 PCFA(1) = P00 * P10 + P10 * P00;              %FROM STATE A2 TO STATE A1        
        %+ P11 * P01 + P01 * P11    %+ P01*P00+P00*P01...
        %+ P00 * P11 + P11 * P00; 
%probability: 2*(1-eps)^3*eps+2*(1-eps)*eps^3+2*(1-eps)^2*eps^2
eps_set_CF_PCFA(1) = 2*(1-eps_set).^3.*eps_set + 2*(1-eps_set).*eps_set.^3 + 2*(1-eps_set).^2.*eps_set.^2; 
 
%Px0C = P00C+P10C; %Px0*Px0;
%Px1C = P01C+P11C; %Px1*Px1;
%}
         
       
RTT = k+M-1;
d = T-RTT;


qsetUncodedARQ = eps_set*r./(1-eps_set);     %This is the GE channel parameter when the system transits into State A1

eps_set_CF = sqrt(eps_set.^4 + 2*eps_set.^3.*(1-eps_set));
qset = eps_set_CF*r./(1-eps_set_CF);         %since eps_B q/(q+r) = eps_CF, and eps_B = 1
                          
eps_set_CF_f = eps_set.^4 + 2*eps_set.^3.*(1-eps_set) + 3*eps_set.^2.*(1-eps_set).^2; %PCF1X(1)=PCF11(1)+PCF10(1)
qf_set = eps_set_CF_f*r./(1-eps_set_CF_f);

eps_set_CF_r = eps_set_CF; %eps_set.^4 + 2*eps_set.^3.*(1-eps_set) + eps_set.^2.*(1-eps_set).^2 + 2*eps_set.*(1-eps_set).^3;  %PCFX1(1)=PCF11(1)+PCF01(1) 
             
qr_set = eps_set_CF_r*r./(1-eps_set_CF_r);
     
%qr_set = qsetUncodedARQ;
 
%{
figure
plot(eps_set_CF)
hold on
plot(eps_set_CF_f,'r')
hold on
plot(eps_set_CF_r,'c')
%}

TransmissionTime = zeros(eps_count,packetno-1);

Delay = zeros(eps_count,packetno-1);

%{
count = 0;
for q = qset
    count = count + 1
    
    qf = qf_set(count);
    qr = qr_set(count);
    
    qUncodedARQ = qsetUncodedARQ(count);
    
    for r_no = 1:total_r_no
        %r_no
        
        %Evolution of the joint GE channel for total_no time slots (Single Realization)
        [ChannelStatef,ChannelStater] = GEchannelCFARQ(qf,qr,r,slotno);        
        %Evolution of GE channel for total_no time slots (Single Realization)
        [ChannelStatef2,ChannelStater2,~] = GEchannelUncodedARQ(qUncodedARQ,r,slotno);
        ChannelStater2 = ChannelStater;
        
        total_no = 10^5;         

        for p = 1:packetno-1 %Each consecutive 2 packets can be considered independently of the others            
            %p
            
            t_tau = M; %time slot for throughput flow graph (need to transmit 2 packets at least once)
            
            t_D = p; %time slot for delay flow graph=(initial) transmission time slot for the p th and p+1 th packets

            total_no = total_no + p;

            %STATES I1 and I2
            ACKreceived = 0;  

            %STATE A2
            t_D = t_D+RTT-1;  
            
            
            %TRANSITION FROM STATE A2 TO STATE B2 (SELF LOOP AROUND A2)
            while sum(ChannelStatef(t_D:t_D+1)) == 0 && max(t_D,t_tau) < 1/3*total_no
                
                t_D = t_D+M;
                if ChannelStater(t_D) == 1 %NACK received (immediate retransmission of both packets)
                   
                    t_D = t_D+RTT-1; 
                    
                %STATE G2 for the DELAY FLOW GRAPH    
                else %NACK is LOST (wait till timeout (d slots) and retransmit 2 packets) 
                    
                    t_D = t_D+d+RTT-1;
                end                
                t_tau = t_tau + M;
            end
            
            %TRANSITION FROM STATE A2 TO STATE A1
            if  (ChannelStatef(t_D-1)+ChannelStater(t_D) == 2 && ChannelStatef(t_D)+ChannelStater(t_D) <= 1) || ... %sum(ChannelStatef(t_D-1:t_D)) == 1 
                (ChannelStatef(t_D)+ChannelStater(t_D) == 2 && ChannelStatef(t_D-1)+ChannelStater(t_D) <= 1) || ...
                (ChannelStatef(t_D-1)+ChannelStatef(t_D) == 1 && ChannelStater(t_D) == 0) 
                         
             %First 2 lines: P00 * P1X + P1X * P00;                     
             %3rd line: P11 * P01 + P01 * P11    
                   
            
                %AT STATE A1, NUMBER OF DoFs ACKNOWLEDGED BY RECEIVER IS N=1    
                %IN THIS CASE, THE REST OF THE PROCESS IS SIMILAR TO UNCODED ARQ
                t_tau = t_tau + 1;      
                t_D = t_D+1;
                
                %TRANSITION FROM STATE A1 TO STATE B1 (SELF LOOP AROUND A1)                 
                while ChannelStatef2(t_D) == 0 && max(t_D,t_tau) < 1/3*total_no %Retransmit packet while forward channel is BAD
                        
                    t_D = t_D+1;
                    if ChannelStater2(t_D) == 1 %NACK received (immediate retransmission)

                        t_D = t_D+k-1;

                    %STATE G1 for the DELAY FLOW GRAPH    
                    else %if max(ChannelStater(t_D+1:t_D+d)) == 0  %NACK is LOST (wait till timeout (d slots) and retransmit packet) 

                        t_D = t_D+d+1+k-1;
                    end    
                    t_tau = t_tau + 1; 
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                t_D_for_tau = t_D;
                %DELAY GRAPH
                %TRANSITION FROM STATE A1 TO STATE 0 

                t_D = t_D+1;
                %TRANSITION FROM STATE A1 TO STATE C1 
                while ChannelStater2(t_D) == 0 && max(t_D,t_tau) < 2/3*total_no  %ACK is LOST                  
                    t_D = t_D+1;                   
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                %THROUGHPUT GRAPH
                t_tau = t_tau + 1; %Need to retransmit 1 of the packets at least once
                while ACKreceived == 0 && max(t_D_for_tau,t_tau) < total_no  %packet will be retransmitted till the ACK is successfully received.

                    %TRANSITION FROM STATE A1 TO STATE C1 THEN TO STATE 0
                    if max(ChannelStater2(t_tau:t_tau+d-1)) == 0 %ACK is LOST and the timer expires   
                        t_tau = t_tau + 1; %not sure about this?
                        while max(ChannelStater2(t_tau:t_tau+T-1)) == 0 && max(t_D_for_tau,t_tau) < total_no 
                            t_tau = t_tau + 1;  %wait till timeout (d+1=T-k slots) and retransmit packet
                        end
                        ACKreceived = 1;


                    %TRANSITION FROM STATE A1 TO STATE 0 
                    elseif max(ChannelStater2(t_tau:t_tau+d-1)) == 1 && max(t_D_for_tau,t_tau) < total_no %ACK arrives before timer expires

                        ACKreceived = 1;

                    end
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %TRANSITION FROM STATE A2 TO STATES 0 OR C2
            else %sum(ChannelStatef(t_D:t_D+1)) == 2 %Forward channel was good for M=2 packets
            
                t_D_for_tau = t_D;
                %DELAY GRAPH
                %TRANSITION FROM STATE A2 TO STATE 0 

                t_D = t_D+1;
                %TRANSITION FROM STATE A2 TO STATE C2 
                while ChannelStater(t_D) == 0 && max(t_D,t_tau) < 2/3*total_no  %ACK is LOST                  
                    t_D = t_D+1;                   
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                %THROUGHPUT GRAPH
                while ACKreceived == 0 && max(t_D_for_tau,t_tau) < total_no  %packet will be retransmitted till the ACK is successfully received.

                    %TRANSITION FROM STATE A2 TO STATE C2 THEN TO STATE 0
                    if max(ChannelStater(t_tau:t_tau+d)) == 0 %ACK is LOST and the timer expires   
                        t_tau = t_tau + M; %not sure about this?
                        while max(ChannelStater(t_tau:t_tau+T-1)) == 0 && max(t_D_for_tau,t_tau) < total_no 
                            t_tau = t_tau + M;  %wait till timeout (d=T-k-1 slots) and retransmit both packets
                        end
                        ACKreceived = 1;


                    %TRANSITION FROM STATE A2 TO STATE 0 
                    elseif max(ChannelStater(t_tau:t_tau+d)) == 1 && max(t_D_for_tau,t_tau) < total_no %ACK arrives before timer expires

                        ACKreceived = 1;

                    end
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                           
           
            end
                                             
            DelayperRealization = t_D-p;
            TransmissionTimeperRealization = t_tau;
            
            Delay(count,p) = Delay(count,p)+DelayperRealization;
            TransmissionTime(count,p) = TransmissionTime(count,p)+TransmissionTimeperRealization;
        end
    end
            
    Delay(count,:) = Delay(count,:)/total_r_no;
    TransmissionTime(count,:) = TransmissionTime(count,:)/total_r_no;
end
%}
 

count = 0;
for q = qf_set
    count = count + 1
    
    qf = qf_set(count);
    qr = qr_set(count);
    
    qUncodedARQ = qsetUncodedARQ(count);
    
    for r_no = 1:total_r_no
        %r_no
        
        %Evolution of the joint GE channel for total_no time slots (Single Realization)
        [ChannelStatef,ChannelStater] = GEchannelCFARQ(qf,qr,r,slotno);        
        %Evolution of GE channel for total_no time slots (Single Realization)
        [ChannelStatef2,ChannelStater2,~] = GEchannelUncodedARQ(qUncodedARQ,r,slotno);
        %ChannelStater2 = ChannelStater;
        
        total_no = 10^5;         

        for p = 1:packetno-1 %Each consecutive 2 packets can be considered independently of the others            
            %p
            
            t_tau = M; %time slot for throughput flow graph (need to transmit 2 packets at least once)
            
            %t_D = p; %time slot for delay flow graph=(initial) transmission time slot for the p th and p+1 th packets
            t_D = randi(packetno,1);
            t_D0 = t_D;
            
            total_no = total_no + t_D0;%p;

            %STATES I1 and I2
            ACKreceived = 0;  

            %STATE A2
            t_D = t_D+RTT-1;  
            
            
            %TRANSITION FROM STATE A2 TO STATE B2 (SELF LOOP AROUND A2)
            while sum(ChannelStatef(t_D:t_D+1)) == 0 && max(t_D,t_tau) < 1/3*total_no
                
                t_D = t_D+M-1;
                if ChannelStater(t_D) == 1 %NACK received (immediate retransmission of both packets)
                   
                    t_D = t_D+RTT-1; 
                    
                %STATE G2 for the DELAY FLOW GRAPH    
                else %NACK is LOST (wait till timeout (d slots) and retransmit 2 packets) 
                    
                    t_D = t_D+d+RTT-1;
                end                
                t_tau = t_tau + M;
            end
            
            %TRANSITION FROM STATE A2 TO STATE A1
            if  (ChannelStatef(t_D-1)+ChannelStater(t_D) == 2 && ChannelStatef(t_D) == 0) || ... %sum(ChannelStatef(t_D-1:t_D)) == 1 
                (ChannelStatef(t_D)  +ChannelStater(t_D) == 2 && ChannelStatef(t_D-1) == 0) || ...
                (ChannelStatef(t_D-1)+ChannelStatef(t_D) == 1 && ChannelStater(t_D) == 0)                        
                %First 2 lines: P00 * P1X + P1X * P00                    
                % and 3rd line: P11 * P01 + P01 * P11    
                   
            
                %AT STATE A1, NUMBER OF DoFs ACKNOWLEDGED BY RECEIVER IS N=1    
                %IN THIS CASE, THE REST OF THE PROCESS IS SIMILAR TO UNCODED ARQ
                t_tau = t_tau + 1;      
                %t_D = t_D+1;
                
                %TRANSITION FROM STATE A1 TO STATE B1 (SELF LOOP AROUND A1)                 
                while ChannelStatef2(t_D) == 0 && max(t_D,t_tau) < 1/3*total_no %Retransmit packet while forward channel is BAD
                        
                    t_D = t_D+1;
                    if ChannelStater2(t_D) == 1 %NACK received (immediate retransmission)

                        t_D = t_D+k-1;

                    %STATE G1 for the DELAY FLOW GRAPH    
                    else %if max(ChannelStater(t_D+1:t_D+d)) == 0  %NACK is LOST (wait till timeout (d slots) and retransmit packet) 

                        t_D = t_D+d+1+k-1;
                    end    
                    t_tau = t_tau + 1; 
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                t_D_for_tau = t_D;
                %DELAY GRAPH
                %TRANSITION FROM STATE A1 TO STATE 0 
                if ChannelStater2(t_D) == 1
                    t_D = t_D+1;
                else
                    %TRANSITION FROM STATE A1 TO STATE C1 
                    while ChannelStater2(t_D) == 0 && max(t_D,t_tau) < 2/3*total_no  %ACK is LOST                  
                        t_D = t_D+1;                   
                    end 
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                %THROUGHPUT GRAPH
                t_tau = t_tau + 1; %Need to retransmit 1 of the packets at least once
                while ACKreceived == 0 && max(t_D_for_tau,t_tau) < total_no  %packet will be retransmitted till the ACK is successfully received.

                    %TRANSITION FROM STATE A1 TO STATE C1 THEN TO STATE 0
                    if max(ChannelStater2(t_tau:t_tau+d-1)) == 0 %ACK is LOST and the timer expires   
                        t_tau = t_tau + 1; %not sure about this?
                        while max(ChannelStater2(t_tau:t_tau+T-1)) == 0 && max(t_D_for_tau,t_tau) < total_no 
                            t_tau = t_tau + 1;  %wait till timeout (d+1=T-k slots) and retransmit packet
                        end
                        ACKreceived = 1;


                    %TRANSITION FROM STATE A1 TO STATE 0 
                    elseif max(ChannelStater2(t_tau:t_tau+d-1)) == 1 && max(t_D_for_tau,t_tau) < total_no %ACK arrives before timer expires

                        ACKreceived = 1;

                    end
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %TRANSITION FROM STATE A2 TO STATES 0 OR C2
            else %sum(ChannelStatef(t_D:t_D+1)) == 2 %Forward channel was good for M=2 packets
            
                t_D_for_tau = t_D;
                %DELAY GRAPH
                %TRANSITION FROM STATE A2 TO STATE 0 
                if ChannelStater(t_D) == 1
                    t_D = t_D+1;
                else
                %TRANSITION FROM STATE A2 TO STATE C2 
                    while ChannelStater(t_D) == 0 && max(t_D,t_tau) < 2/3*total_no  %ACK is LOST                  
                        t_D = t_D+1;                   
                    end
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                %THROUGHPUT GRAPH
                while ACKreceived == 0 && max(t_D_for_tau,t_tau) < total_no  %packet will be retransmitted till the ACK is successfully received.

                    %TRANSITION FROM STATE A2 TO STATE C2 THEN TO STATE 0
                    if max(ChannelStater(t_tau:t_tau+d)) == 0 %ACK is LOST and the timer expires   
                        t_tau = t_tau + M; %not sure about this?
                        while max(ChannelStater(t_tau:t_tau+T-1)) == 0 && max(t_D_for_tau,t_tau) < total_no 
                            t_tau = t_tau + M;  %wait till timeout (d=T-k-1 slots) and retransmit both packets
                        end
                        ACKreceived = 1;


                    %TRANSITION FROM STATE A2 TO STATE 0 
                    elseif max(ChannelStater(t_tau:t_tau+d)) == 1 && max(t_D_for_tau,t_tau) < total_no %ACK arrives before timer expires

                        ACKreceived = 1;

                    end
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                           
           
            end
                                             
            DelayperRealization = t_D-t_D0;%p;
            TransmissionTimeperRealization = t_tau;
            
            Delay(count,p) = Delay(count,p)+DelayperRealization;
            TransmissionTime(count,p) = TransmissionTime(count,p)+TransmissionTimeperRealization;
        end
    end
            
    Delay(count,:) = Delay(count,:)/total_r_no;
    TransmissionTime(count,:) = TransmissionTime(count,:)/total_r_no;
end

MeanDelay = zeros(1,eps_count);
ThroughputInv = zeros(1,eps_count);
for i = 1:eps_count
    for j = 1:packetno-1                
        MeanDelay(i) = MeanDelay(i) + Delay(i,j);
        ThroughputInv(i) = ThroughputInv(i) + TransmissionTime(i,j);
    end
    MeanDelay(i) = MeanDelay(i)/packetno;
    ThroughputInv(i) = ThroughputInv(i)/packetno;
end
Throughput = M./ThroughputInv;
    
t = 1;
T_count = 1;
clear str;  str = cell(1,T_count);
str{t} = ['CF ARQ, T=' num2str(T)]; 

figure    
plot(eps_set,MeanDelay,'k','linewidth',2,...
        'color',colorStyles{1+rem(t,numel(colorStyles))},'linestyle',lineStyles{1+rem(t,numel(lineStyles))},...
        'marker',markerStyles{1+rem(t,numel(markerStyles))},'markersize',10); hold on;  
xlab = 'Erasure rate, $\epsilon$'; 
ylab = 'Average delay, $\bar{D}$';
box on;     set(gca,'FontSize',20) 
xlhand = get(gca,'xlabel'); xlabel(xlab,'Interpreter','latex'); set(xlhand,'fontsize',20) 
ylhand = get(gca,'ylabel'); ylabel(ylab,'Interpreter','latex'); set(ylhand,'fontsize',20)
legend(str,'FontSize',20,'location','NorthEast');
xaxis = 0.3;    yaxis = 5;
text(xaxis,yaxis,['M=' num2str(M) ', N=' num2str(N) ', r=' num2str(r) ', k=' num2str(k)],'fontsize',20);
   
figure
plot(eps_set,Throughput,'k','linewidth',2,...
        'color',colorStyles{1+rem(t,numel(colorStyles))},'linestyle',lineStyles{1+rem(t,numel(lineStyles))},...
        'marker',markerStyles{1+rem(t,numel(markerStyles))},'markersize',10); hold on;    
xlab = 'Erasure rate, $\epsilon$'; 
ylab = 'Throughput, $\eta$';
box on;     set(gca,'FontSize',20) 
xlhand = get(gca,'xlabel'); xlabel(xlab,'Interpreter','latex'); set(xlhand,'fontsize',20) 
ylhand = get(gca,'ylabel'); ylabel(ylab,'Interpreter','latex'); set(ylhand,'fontsize',20)
legend(str,'FontSize',20,'location','NorthEast');
xaxis = 0.3;    yaxis = 0.35;
text(xaxis,yaxis,['M=' num2str(M) ', N=' num2str(N) ', r=' num2str(r) ', k=' num2str(k)],'fontsize',20);


%% 4-Coded ARQ scheme with no HARQ combining 

%M=2 packets and N=2

RTT = k+M-1;
d = T-RTT;

eps_set_CF_f = eps_set.^4+2*eps_set.^3.*(1-eps_set) + 3*eps_set.^2.*(1-eps_set).^2;
qf_set = eps_set_CF_f*r./(1-eps_set_CF_f);
%eps_set_CF_r = eps_set.^4+2*eps_set.^3.*(1-eps_set) + (1-eps_set).^2.*eps_set.^2 + 2*(1-eps_set).^3.*eps_set;
%qr_set = eps_set_CF_r*r./(1-eps_set_CF_r);
qr_set = eps_set*r./(1-eps_set);

eps_set_CF = sqrt(eps_set.^4+2*eps_set.^3.*(1-eps_set));
qset = eps_set_CF*r./(1-eps_set_CF);  %since eps_B q/(q+r) = eps_CF, and eps_B = 1
qsetUncodedARQ = eps_set*r./(1-eps_set);  %This is the GE channel parameter when the system transits into State A1

TransmissionTime = zeros(eps_count,packetno-1);

Delay = zeros(eps_count,packetno-1);

%{
count = 0;
for q = qset
    count = count + 1
    
    qf = qf_set(count);
    qr = qr_set(count);
    
    qUncodedARQ = qsetUncodedARQ(count);
    
    for r_no = 1:total_r_no
              
        %Evolution of the joint GE channel for total_no time slots (Single Realization)
        [ChannelStatef,ChannelStater] = GEchannelCFARQ(qf,qr,r,slotno);        
        %Evolution of GE channel for total_no time slots (Single Realization)
        [ChannelStatef2,ChannelStater2,~] = GEchannelUncodedARQ(qUncodedARQ,r,slotno);
        ChannelStater2 = ChannelStater;
        
        
        total_no = 10^5;
        
        for p = 1:packetno-1 %Each consecutive 2 packets can be considered independently of the others            
            %p
            
            t_tau = M; %time slot for throughput flow graph (need to transmit 2 packets at least once)
            
            t_D = p; %time slot for delay flow graph=(initial) transmission time slot for the p th and p+1 th packets

            total_no = total_no + p;

            %STATE I1 and I2
            ACKreceived = 0;  

            %STATE A2
            t_D = t_D+RTT-1;  
            
            %TRANSITION FROM STATE A2 TO STATE B2 (SELF LOOP AROUND A2)
            while sum(ChannelStatef(t_D:t_D+1)) == 0 && max(t_D,t_tau) < 1/3*total_no
                
                t_D = t_D+M;
                if ChannelStater(t_D) == 1 %NACK received (immediate retransmission of both packets)
                   
                    t_D = t_D+RTT-1; 
                    
                %STATE G2 for the DELAY FLOW GRAPH    
                else %NACK is LOST (wait till timeout (d slots) and retransmit 2 packets) 
                    
                    t_D = t_D+d+RTT-1;
                end                
                t_tau = t_tau + M;
            end
            
            %TRANSITION FROM STATE A2 TO STATES G3 AND A3
            if sum(ChannelStatef(t_D:t_D+1)) == 1 && ChannelStater(t_D) == 0
                
                %TRANSITION FROM STATE G3 TO STATE B3 THEN TO STATE A3
                t_D = t_D+M;
                
                %NACK is LOST (wait till timeout (d slots) and retransmit 2 packets)
                while ChannelStatef2(t_D) == 0 && ChannelStater2(t_D) == 0 && max(t_D,t_tau) < 1/3*total_no
                
                    t_D = t_D+M;                   

                    %STATE G3 for the DELAY FLOW GRAPH    
                    t_D = t_D+M*(d+RTT-1);
                                    
                    t_tau = t_tau + M; %both packets are retransmitted here
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
                              
                t_D_init = t_D;
                %TRANSITION FROM STATE A3 TO STATE A1
                if ChannelStatef2(t_D) == 0 && ChannelStater2(t_D) == 1
                
                    %AT STATE A1, NUMBER OF DoFs ACKNOWLEDGED BY RECEIVER IS N=1    
                    %IN THIS CASE, THE REST OF THE PROCESS IS SIMILAR TO UNCODED ARQ
                    %t_tau = t_tau + 1;      
                    t_D = t_D+1;

                    %TRANSITION FROM STATE A1 TO STATE B1 (SELF LOOP AROUND A1)                 
                    while ChannelStatef2(t_D) == 0 && max(t_D,t_tau) < 1/3*total_no %Retransmit packet while forward channel is BAD

                        t_D = t_D+1;
                        if ChannelStater2(t_D) == 1 %NACK received (immediate retransmission)

                            t_D = t_D+k-1;

                        %STATE G1 for the DELAY FLOW GRAPH    
                        else %if max(ChannelStater(t_D+1:t_D+d)) == 0  %NACK is LOST (wait till timeout (d slots) and retransmit packet) 

                            t_D = t_D+d+1+k-1;
                        end    
                        t_tau = t_tau + 1; 
                    end
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    t_D_for_tau = t_D;
                    %DELAY GRAPH
                    %TRANSITION FROM STATE A1 TO STATE 0 

                    t_D = t_D+1;
                    %TRANSITION FROM STATE A1 TO STATE C1 
                    while ChannelStater2(t_D) == 0 && max(t_D,t_tau) < 2/3*total_no  %ACK is LOST                  
                        t_D = t_D+1;                   
                    end
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    %THROUGHPUT GRAPH
                    t_tau = t_tau + 1; %Need to retransmit 1 of the packets at least once
                    while ACKreceived == 0 && max(t_D_for_tau,t_tau) < total_no  %packet will be retransmitted till the ACK is successfully received.

                        %TRANSITION FROM STATE A1 TO STATE C1 THEN TO STATE 0
                        if max(ChannelStater2(t_tau:t_tau+d-1)) == 0 %ACK is LOST and the timer expires   
                            t_tau = t_tau + 1; %not sure about this?
                            while max(ChannelStater2(t_tau:t_tau+T-1)) == 0 && max(t_D_for_tau,t_tau) < total_no 
                                t_tau = t_tau + 1;  %wait till timeout (d+1=T-k slots) and retransmit packet
                            end
                            ACKreceived = 1;


                        %TRANSITION FROM STATE A1 TO STATE 0 
                        elseif max(ChannelStater2(t_tau:t_tau+d-1)) == 1 && max(t_D_for_tau,t_tau) < total_no %ACK arrives before timer expires

                            ACKreceived = 1;

                        end
                    end
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                             
                    
                %TRANSITION FROM STATE A3 TO STATES C3 AND 0 (SAME AS THE TRANSITION FROM STATE A1 TO C1 AND 0)
                else
                    
                    t_D_for_tau = t_D;
                    %DELAY GRAPH
                    %TRANSITION FROM STATE A1 TO STATE 0 

                    t_D = t_D+1;
                    %TRANSITION FROM STATE A1 TO STATE C1 
                    while ChannelStater2(t_D) == 0 && max(t_D,t_tau) < 2/3*total_no  %ACK is LOST                  
                        t_D = t_D+1;                   
                    end
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    %THROUGHPUT GRAPH
                    t_tau = t_tau + 1; %Need to retransmit 1 of the packets at least once
                    while ACKreceived == 0 && max(t_D_for_tau,t_tau) < total_no  %packet will be retransmitted till the ACK is successfully received.

                        %TRANSITION FROM STATE A1 TO STATE C1 THEN TO STATE 0
                        if max(ChannelStater2(t_tau:t_tau+d-1)) == 0 %ACK is LOST and the timer expires   
                            %t_tau = t_tau + 1; %not sure about this?
                            while max(ChannelStater2(t_tau:t_tau+T-1)) == 0 && max(t_D_for_tau,t_tau) < total_no 
                                t_tau = t_tau + 1;  %wait till timeout (d+1=T-k slots) and retransmit packet
                            end
                            ACKreceived = 1;


                        %TRANSITION FROM STATE A1 TO STATE 0 
                        elseif max(ChannelStater2(t_tau:t_tau+d-1)) == 1 && max(t_D_for_tau,t_tau) < total_no %ACK arrives before timer expires

                            ACKreceived = 1;

                        end
                    end
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
                
                end
                gap = t_D - t_D_init;
                t_D = t_D_init+2*gap;
            
            %TRANSITION FROM STATE A2 TO STATE A1
            elseif sum(ChannelStatef(t_D:t_D+1)) == 1 && ChannelStater(t_D) == 1
                
                %AT STATE A1, NUMBER OF DoFs ACKNOWLEDGED BY RECEIVER IS N=1    
                %IN THIS CASE, THE REST OF THE PROCESS IS SIMILAR TO UNCODED ARQ
                %t_tau = t_tau + 1;      
                t_D = t_D+1;
                
                %TRANSITION FROM STATE A1 TO STATE B1 (SELF LOOP AROUND A1)                 
                while ChannelStatef2(t_D) == 0 && max(t_D,t_tau) < 1/3*total_no %Retransmit packet while forward channel is BAD
                        
                    t_D = t_D+1;
                    if ChannelStater2(t_D) == 1 %NACK received (immediate retransmission)

                        t_D = t_D+k-1;

                    %STATE G1 for the DELAY FLOW GRAPH    
                    else %if max(ChannelStater(t_D+1:t_D+d)) == 0  %NACK is LOST (wait till timeout (d slots) and retransmit packet) 

                        t_D = t_D+d+1+k-1;
                    end    
                    t_tau = t_tau + 1; 
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                t_D_for_tau = t_D;
                %DELAY GRAPH
                %TRANSITION FROM STATE A1 TO STATE 0 

                t_D = t_D+1;
                %TRANSITION FROM STATE A1 TO STATE C1 
                while ChannelStater2(t_D) == 0 && max(t_D,t_tau) < 2/3*total_no  %ACK is LOST                  
                    t_D = t_D+1;                   
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                %THROUGHPUT GRAPH
                t_tau = t_tau + 1; %Need to retransmit 1 of the packets at least once
                while ACKreceived == 0 && max(t_D_for_tau,t_tau) < total_no  %packet will be retransmitted till the ACK is successfully received.

                    %TRANSITION FROM STATE A1 TO STATE C1 THEN TO STATE 0
                    if max(ChannelStater2(t_tau:t_tau+d-1)) == 0 %ACK is LOST and the timer expires   
                        t_tau = t_tau + 1; %not sure about this?
                        while max(ChannelStater2(t_tau:t_tau+T-1)) == 0 && max(t_D_for_tau,t_tau) < total_no 
                            t_tau = t_tau + 1;  %wait till timeout (d+1=T-k slots) and retransmit packet
                        end
                        ACKreceived = 1;


                    %TRANSITION FROM STATE A1 TO STATE 0 
                    elseif max(ChannelStater2(t_tau:t_tau+d-1)) == 1 && max(t_D_for_tau,t_tau) < total_no %ACK arrives before timer expires

                        ACKreceived = 1;

                    end
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%             

            %TRANSITION FROM STATE A2 TO STATES 0 OR C2
            else
                
                t_D_for_tau = t_D;
                %DELAY GRAPH
                %TRANSITION FROM STATE A2 TO STATE 0 

                t_D = t_D+1;
                %TRANSITION FROM STATE A2 TO STATE C2 
                while ChannelStater(t_D) == 0 && max(t_D,t_tau) < 2/3*total_no  %ACK is LOST                  
                    t_D = t_D+1;                   
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                %THROUGHPUT GRAPH
                while ACKreceived == 0 && max(t_D_for_tau,t_tau) < total_no  %packet will be retransmitted till the ACK is successfully received.

                    %TRANSITION FROM STATE A2 TO STATE C2 THEN TO STATE 0
                    if max(ChannelStater(t_tau:t_tau+d)) == 0 %ACK is LOST and the timer expires   
                        t_tau = t_tau + M; %not sure about this?
                        while max(ChannelStater(t_tau:t_tau+T-1)) == 0 && max(t_D_for_tau,t_tau) < total_no 
                            t_tau = t_tau + M;  %wait till timeout (d=T-k-1 slots) and retransmit both packets
                        end
                        ACKreceived = 1;


                    %TRANSITION FROM STATE A2 TO STATE 0 
                    elseif max(ChannelStater(t_tau:t_tau+d)) == 1 && max(t_D_for_tau,t_tau) < total_no %ACK arrives before timer expires

                        ACKreceived = 1;

                    end
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                          
                
            end
            
            
            DelayperRealization = t_D-p;
            TransmissionTimeperRealization = t_tau;
            
            Delay(count,p) = Delay(count,p)+DelayperRealization;
            TransmissionTime(count,p) = TransmissionTime(count,p)+TransmissionTimeperRealization;
        end
    
    
    end
    Delay(count,:) = Delay(count,:)/total_r_no;
    TransmissionTime(count,:) = TransmissionTime(count,:)/total_r_no;    
end
%}

count = 0;
for q = qset
    count = count + 1
    
    qf = qf_set(count);
    qr = qr_set(count);
    
    qUncodedARQ = qsetUncodedARQ(count);
    
    for r_no = 1:total_r_no
              
        %Evolution of the joint GE channel for total_no time slots (Single Realization)
        [ChannelStatef,ChannelStater] = GEchannelCFARQ(qf,qr,r,slotno);        
        %Evolution of GE channel for total_no time slots (Single Realization)
        [ChannelStatef2,ChannelStater2,~] = GEchannelUncodedARQ(qUncodedARQ,r,slotno);
        ChannelStater2 = ChannelStater;
        
        
        total_no = 10^5;
        
        for p = 1:packetno-1 %Each consecutive 2 packets can be considered independently of the others            
            %p
            
            t_tau = M; %time slot for throughput flow graph (need to transmit 2 packets at least once)
            
            t_D = p; %time slot for delay flow graph=(initial) transmission time slot for the p th and p+1 th packets

            total_no = total_no + p;

            %STATE I1 and I2
            ACKreceived = 0;  

            %STATE A2
            t_D = t_D+RTT-1;  
            
            %TRANSITION FROM STATE A2 TO STATE B2 (SELF LOOP AROUND A2)
            while sum(ChannelStatef(t_D:t_D+1)) == 0 && max(t_D,t_tau) < 1/3*total_no
                
                t_D = t_D+M;
                if ChannelStater(t_D) == 1 %NACK received (immediate retransmission of both packets)
                   
                    t_D = t_D+RTT-1; 
                    
                %STATE G2 for the DELAY FLOW GRAPH    
                else %NACK is LOST (wait till timeout (d slots) and retransmit 2 packets) 
                    
                    t_D = t_D+d+RTT-1;
                end                
                t_tau = t_tau + M;
            end
            
            %TRANSITION FROM STATE A2 TO STATES G3 AND A3
            if sum(ChannelStatef(t_D:t_D+1)) == 1 && ChannelStater(t_D) == 0
                
                %TRANSITION FROM STATE G3 TO STATE B3 THEN TO STATE A3
                t_D = t_D+M;
                
                %NACK is LOST (wait till timeout (d slots) and retransmit 2 packets)
                while ChannelStatef2(t_D) == 0 && ChannelStater2(t_D) == 0 && max(t_D,t_tau) < 1/3*total_no
                
                    t_D = t_D+M;                   

                    %STATE G3 for the DELAY FLOW GRAPH    
                    t_D = t_D+M*(d+RTT-1);
                                    
                    t_tau = t_tau + 1; %one packet is retransmitted here
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
                              
                t_D_init = t_D;
                %TRANSITION FROM STATE A3 TO STATE A1
                if ChannelStatef2(t_D) == 0 && ChannelStater2(t_D) == 1
                
                    %AT STATE A1, NUMBER OF DoFs ACKNOWLEDGED BY RECEIVER IS N=1    
                    %IN THIS CASE, THE REST OF THE PROCESS IS SIMILAR TO UNCODED ARQ
                    %t_tau = t_tau + 1;      
                    t_D = t_D+1;

                    %TRANSITION FROM STATE A1 TO STATE B1 (SELF LOOP AROUND A1)                 
                    while ChannelStatef2(t_D) == 0 && max(t_D,t_tau) < 1/3*total_no %Retransmit packet while forward channel is BAD

                        t_D = t_D+1;
                        if ChannelStater2(t_D) == 1 %NACK received (immediate retransmission)

                            t_D = t_D+k-1;

                        %STATE G1 for the DELAY FLOW GRAPH    
                        else %if max(ChannelStater(t_D+1:t_D+d)) == 0  %NACK is LOST (wait till timeout (d slots) and retransmit packet) 

                            t_D = t_D+d+1+k-1;
                        end    
                        t_tau = t_tau + 1; 
                    end
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    t_D_for_tau = t_D;
                    %DELAY GRAPH
                    %TRANSITION FROM STATE A1 TO STATE 0 

                    t_D = t_D+1;
                    %TRANSITION FROM STATE A1 TO STATE C1 
                    while ChannelStater2(t_D) == 0 && max(t_D,t_tau) < 2/3*total_no  %ACK is LOST                  
                        t_D = t_D+1;                   
                    end
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    %THROUGHPUT GRAPH
                    while ACKreceived == 0 && max(t_D_for_tau,t_tau) < total_no  %packet will be retransmitted till the ACK is successfully received.

                        %TRANSITION FROM STATE A1 TO STATE C1 THEN TO STATE 0
                        if max(ChannelStater2(t_tau:t_tau+d-1)) == 0 %ACK is LOST and the timer expires   
                            while max(ChannelStater2(t_tau:t_tau+T-1)) == 0 && max(t_D_for_tau,t_tau) < total_no 
                                t_tau = t_tau + 1;  %wait till timeout (d+1=T-k slots) and retransmit packet
                            end
                            ACKreceived = 1;


                        %TRANSITION FROM STATE A1 TO STATE 0 
                        elseif max(ChannelStater2(t_tau:t_tau+d-1)) == 1 && max(t_D_for_tau,t_tau) < total_no %ACK arrives before timer expires

                            ACKreceived = 1;

                        end
                    end
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                             
                    
                %TRANSITION FROM STATE A3 TO STATES C3 AND 0 (SAME AS THE TRANSITION FROM STATE A1 TO C1 AND 0)
                else
                    
                    t_D_for_tau = t_D;
                    %DELAY GRAPH
                    %TRANSITION FROM STATE A1 TO STATE 0 

                    t_D = t_D+1;
                    %TRANSITION FROM STATE A1 TO STATE C1 
                    while ChannelStater2(t_D) == 0 && max(t_D,t_tau) < 2/3*total_no  %ACK is LOST                  
                        t_D = t_D+1;                   
                    end
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    %THROUGHPUT GRAPH
                    while ACKreceived == 0 && max(t_D_for_tau,t_tau) < total_no  %packet will be retransmitted till the ACK is successfully received.

                        %TRANSITION FROM STATE A1 TO STATE C1 THEN TO STATE 0
                        if max(ChannelStater2(t_tau:t_tau+d-1)) == 0 %ACK is LOST and the timer expires   
                            while max(ChannelStater2(t_tau:t_tau+T-1)) == 0 && max(t_D_for_tau,t_tau) < total_no 
                                t_tau = t_tau + 1;  %wait till timeout (d+1=T-k slots) and retransmit packet
                            end
                            ACKreceived = 1;


                        %TRANSITION FROM STATE A1 TO STATE 0 
                        elseif max(ChannelStater2(t_tau:t_tau+d-1)) == 1 && max(t_D_for_tau,t_tau) < total_no %ACK arrives before timer expires

                            ACKreceived = 1;

                        end
                    end
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
                
                end
                gap = t_D - t_D_init;
                t_D = t_D_init+2*gap;
            
            %TRANSITION FROM STATE A2 TO STATE A1
            elseif sum(ChannelStatef(t_D:t_D+1)) == 1 && ChannelStater(t_D) == 1
                
                %AT STATE A1, NUMBER OF DoFs ACKNOWLEDGED BY RECEIVER IS N=1    
                %IN THIS CASE, THE REST OF THE PROCESS IS SIMILAR TO UNCODED ARQ
                t_D = t_D+1;
                
                %TRANSITION FROM STATE A1 TO STATE B1 (SELF LOOP AROUND A1)                 
                while ChannelStatef2(t_D) == 0 && max(t_D,t_tau) < 1/3*total_no %Retransmit packet while forward channel is BAD
                        
                    t_D = t_D+1;
                    if ChannelStater2(t_D) == 1 %NACK received (immediate retransmission)

                        t_D = t_D+k-1;

                    %STATE G1 for the DELAY FLOW GRAPH    
                    else %if max(ChannelStater(t_D+1:t_D+d)) == 0  %NACK is LOST (wait till timeout (d slots) and retransmit packet) 

                        t_D = t_D+d+1+k-1;
                    end    
                    t_tau = t_tau + 1; 
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                t_D_for_tau = t_D;
                %DELAY GRAPH
                %TRANSITION FROM STATE A1 TO STATE 0 

                t_D = t_D+1;
                %TRANSITION FROM STATE A1 TO STATE C1 
                while ChannelStater2(t_D) == 0 && max(t_D,t_tau) < 2/3*total_no  %ACK is LOST                  
                    t_D = t_D+1;                   
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                %THROUGHPUT GRAPH
                while ACKreceived == 0 && max(t_D_for_tau,t_tau) < total_no  %packet will be retransmitted till the ACK is successfully received.

                    %TRANSITION FROM STATE A1 TO STATE C1 THEN TO STATE 0
                    if max(ChannelStater2(t_tau:t_tau+d-1)) == 0 %ACK is LOST and the timer expires   
                        while max(ChannelStater2(t_tau:t_tau+T-1)) == 0 && max(t_D_for_tau,t_tau) < total_no 
                            t_tau = t_tau + 1;  %wait till timeout (d+1=T-k slots) and retransmit packet
                        end
                        ACKreceived = 1;


                    %TRANSITION FROM STATE A1 TO STATE 0 
                    elseif max(ChannelStater2(t_tau:t_tau+d-1)) == 1 && max(t_D_for_tau,t_tau) < total_no %ACK arrives before timer expires

                        ACKreceived = 1;

                    end
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%             

            %TRANSITION FROM STATE A2 TO STATES 0 OR C2
            else
                
                t_D_for_tau = t_D;
                %DELAY GRAPH
                %TRANSITION FROM STATE A2 TO STATE 0 

                t_D = t_D+1;
                %TRANSITION FROM STATE A2 TO STATE C2 
                while ChannelStater(t_D) == 0 && max(t_D,t_tau) < 2/3*total_no  %ACK is LOST                  
                    t_D = t_D+1;                   
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                %THROUGHPUT GRAPH
                while ACKreceived == 0 && max(t_D_for_tau,t_tau) < total_no  %packet will be retransmitted till the ACK is successfully received.

                    %TRANSITION FROM STATE A2 TO STATE C2 THEN TO STATE 0
                    if max(ChannelStater(t_tau:t_tau+d)) == 0 %ACK is LOST and the timer expires   
                        while max(ChannelStater(t_tau:t_tau+T-1)) == 0 && max(t_D_for_tau,t_tau) < total_no 
                            t_tau = t_tau + M;  %wait till timeout (d=T-k-1 slots) and retransmit both packets
                        end
                        ACKreceived = 1;


                    %TRANSITION FROM STATE A2 TO STATE 0 
                    elseif max(ChannelStater(t_tau:t_tau+d)) == 1 && max(t_D_for_tau,t_tau) < total_no %ACK arrives before timer expires

                        ACKreceived = 1;

                    end
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                          
                
            end
            
            
            DelayperRealization = t_D-p;
            TransmissionTimeperRealization = t_tau;
            
            Delay(count,p) = Delay(count,p)+DelayperRealization;
            TransmissionTime(count,p) = TransmissionTime(count,p)+TransmissionTimeperRealization;
        end
    
    
    end
    Delay(count,:) = Delay(count,:)/total_r_no;
    TransmissionTime(count,:) = TransmissionTime(count,:)/total_r_no;    
end

MeanDelay = zeros(1,eps_count);
ThroughputInv = zeros(1,eps_count);
for i = 1:eps_count
    for j = 1:packetno-1                
        MeanDelay(i) = MeanDelay(i) + Delay(i,j);
        ThroughputInv(i) = ThroughputInv(i) + TransmissionTime(i,j);
    end
    MeanDelay(i) = MeanDelay(i)/packetno;
    ThroughputInv(i) = ThroughputInv(i)/packetno;
end
Throughput = M./ThroughputInv;

t = 1;
T_count = 1;
clear str;  str = cell(1,T_count);
str{t} = ['C-ARQ, T=' num2str(T)]; 

figure    
plot(eps_set,MeanDelay,'k','linewidth',2,...
        'color',colorStyles{1+rem(t,numel(colorStyles))},'linestyle',lineStyles{1+rem(t,numel(lineStyles))},...
        'marker',markerStyles{1+rem(t,numel(markerStyles))},'markersize',10); hold on;  
xlab = 'Erasure rate, $\epsilon$'; 
ylab = 'Average delay, $\bar{D}$';
box on;     set(gca,'FontSize',20) 
xlhand = get(gca,'xlabel'); xlabel(xlab,'Interpreter','latex'); set(xlhand,'fontsize',20) 
ylhand = get(gca,'ylabel'); ylabel(ylab,'Interpreter','latex'); set(ylhand,'fontsize',20)
legend(str,'FontSize',20,'location','NorthEast');
xaxis = 0.3;    yaxis = 5;
text(xaxis,yaxis,['M=' num2str(M) ', N=' num2str(N) ', r=' num2str(r) ', k=' num2str(k)],'fontsize',20);
   
figure
plot(eps_set,Throughput,'k','linewidth',2,...
        'color',colorStyles{1+rem(t,numel(colorStyles))},'linestyle',lineStyles{1+rem(t,numel(lineStyles))},...
        'marker',markerStyles{1+rem(t,numel(markerStyles))},'markersize',10); hold on;    
xlab = 'Erasure rate, $\epsilon$'; 
ylab = 'Throughput, $\eta$';
box on;     set(gca,'FontSize',20) 
xlhand = get(gca,'xlabel'); xlabel(xlab,'Interpreter','latex'); set(xlhand,'fontsize',20) 
ylhand = get(gca,'ylabel'); ylabel(ylab,'Interpreter','latex'); set(ylhand,'fontsize',20)
legend(str,'FontSize',20,'location','NorthEast');
xaxis = 0.3;    yaxis = 0.35;
text(xaxis,yaxis,['M=' num2str(M) ', N=' num2str(N) ', r=' num2str(r) ', k=' num2str(k)],'fontsize',20);


%% 5-Coded ARQ scheme with HARQ combining
%This is left as future work.



