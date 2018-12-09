function [phiD, PMF_delay, CCDF_delay] = HARQtailDelay(k,K,T,P_kron,pi_I_kron,P00,P01,P10,P11,P0x,Px0,P1x,Px1,r,NACK,epsf,eps,Ztransform)

    syms z
    
    d = T-k;      %time to timeout
    t = T; 

    eps_G0 = 0; %0.07;
    eps_B0 = 1; %0.7;
    if K~=1
        I = eye(2*K,2*K);
        onevec = ones(2*K,1);
        zeromat = zeros(2*K,2*K);
        onemat = ones(2*K,2*K);
    else %K==1
        I = 1;
        onevec = 1;
        zeromat = 0;
        onemat = 1;
    end

    %{
    %ARQ
    GammaGperRho = -log(1-eps_G0);
    GammaBperRho = -log(1-eps_B0);
    
    eps_G_function = @(m) 1-exp(-GammaGperRho/m);
    eps_B_function = @(m) 1-exp(-GammaBperRho/m);
    %}
    
    %%{
    %HARQ
    GammaperRho = 10*eps;%-log(1-eps)*5;    
    eps_G_function = @(m) eps_G0*(1-exp(-GammaperRho./m));
    eps_B_function = @(m) eps_B0*(1-exp(-GammaperRho./m));
    
    %figure; plot(1:100,eps_B_function(1:100))
    %}       
    
    
    prod1 = I;
    D1 = zeromat;
    D2 = zeromat;
    D3 = onemat;
    
    C1 = zeromat;
    C2 = zeromat;
    
    sumtemp = I;
    
    cmax = max(T,10);       % 100;  
    
    
    for c = 0:cmax   

       if c == 0
           
          [P00,P01,P10,P11,P0x,Px0,P1x,Px1,P_kron,pi_kron,pi_I_kron] = TransitionProbabilityMatrices(eps_G_function(1),eps_B_function(1),eps,eps_G_function(1),eps_B_function(1),eps,r,r,NACK,K);

           D1 = Px0;                    %Px0(k+1)
           
           %
           prodc1 = z^k*P10*P_kron^(k-1)-z^T*P11*P_kron^(T-1);
           prodc1derivative = k*z^(k-1)*P10*P_kron^(k-1)-T*z^(T-1)*P11*P_kron^(T-1);
           C1derivative = prodc1derivative;
           prodc2 = I*Px0;           
           C2derivative = zeromat;
           %
           
           
           if d == 1
               sum1 = D1;
           elseif d== 0
               sum1 = zeromat;    
           end
            
       else %c >= 1

           eps_G = eps_G_function(c);
           eps_B = eps_B_function(c);        
           [P00,P01,P10,P11,P0x,Px0,P1x,Px1,P_kron,pi_kron,pi_I_kron] = TransitionProbabilityMatrices(eps_G,eps_B,eps,eps_G,eps_B,eps,r,r,NACK,K);

           prod1 = prod1 * Px1;         %prod_{i=0}^k {Px1(i)}   

           %
           prodc1 = prodc1 * (z^k*P10*P_kron^(k-1)-z^T*P11*P_kron^(T-1));
           prodc1derivative = prodc1derivative * (z^k*P10*P_kron^(k-1)-z^T*P11*P_kron^(T-1)); %this is wrong, check pls!
           
           prodc2 = prodc2 * Px1;  
           [P00,P01,P10,P11,P0x,Px0,P1x,Px1,P_kron,pi_kron,pi_I_kron] = TransitionProbabilityMatrices(eps_G_function(c+1),eps_B_function(c+1),eps,eps_G_function(c+1),eps_B_function(c+1),eps,r,r,NACK,K);
           prodc2 = prodc2 * Px0; 
           %
           
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           if d >= 1
               if c == d 
                    sumtemp = prod1; 
               end            
           end
           
           if c == T              
                %tempmat = I / (I-z*prod1); %inverse of a matrix (that can be close to singular)           
                tempmat = zeromat;
                for j = 0:10
                    tempmat =  tempmat + (z*prod1)^j;
                end
                
                sum2 = sumtemp * tempmat; %This is the 2nd term in throughput analysis 
                sum2derivative = sumtemp * tempmat * prod1 * tempmat;
           end
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

           
           eps_G = eps_G_function(c+1);
           eps_B = eps_B_function(c+1);        
           [P00,P01,P10,P11,P0x,Px0,P1x,Px1,P_kron,pi_kron,pi_I_kron] = TransitionProbabilityMatrices(eps_G,eps_B,eps,eps_G,eps_B,eps,r,r,NACK,K);

           D1 = D1 + prod1 * Px0;       %Px0(k+1)
           D2 = D2 + c * prod1 * Px0; 
           D3 = D3 * prod1 * Px0;
           
           %
           C1 = C1 + prodc1;
           C1derivative = C1derivative + prodc1derivative;
           C2 = C2 + z^c*prodc2;
           C2derivative = C2derivative + c*z^(c-1)*prodc2;
           %
           
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           if d > 1 
               if c == d-1
                   sum1 = D1; %This is the 1st term in throughput analysis
               end
           end
           
           if c == T-1
              sum3 = D1;
           end
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


       end


    end  
 
    
    [P00,P01,P10,P11,P0x,Px0,P1x,Px1,P_kron,pi_kron,pi_I_kron] = TransitionProbabilityMatrices(eps_G0,eps_B0,eps,eps_G0,eps_B0,eps,r,r,NACK,K);
    
    A = @(z) I-z^k*P10*P_kron^(k-1)-z^T*P11*P_kron^(T-1);
    %B = I-P10*P_kron^(k-1)-P11*P_kron^(t-1);
    %C = I-Px1;
    
    
    %b*inv(A) = b/A         %inv(A)*b = A\b
    
    
    % ANALYTICAL EXPRESSION FOR DELAY

    PhiD = @(z) z^(k-1)*P_kron^(k-1)/A(z)*(z*P00+z^2*(P01/(I-z*Px1)) *Px0);
    %PhiD = z^(k-1)*P_kron^(k-1)*C1.*(z*P00+z^2*P01*C2); %updated
    
             
    phiD = @(z) (pi_I_kron * PhiD(z) * onevec)/(pi_I_kron * onevec);
    
    if Ztransform==1
        %Z-transform of x[n]: sum_n=-inf^inf x[n]z^-n
        %If x[n]=P(X=n), and z=z^-1, this gives the probability-generating function (PGF) for X: E[z^X].
        n_min = k;
        n_max = 5*T;
        no = n_max-n_min+1;
        PMF_delay = zeros(1,no);
        for count = 1:no
            count
            n = n_min+count-1;
            PMF_delay(count) = iztrans(phiD(1/z),n);
        end
        figure
        plot(n_min:n_max,PMF_delay,'k','linewidth',2)
        xlab = 'Value of delay, d'; 
        ylab = 'Probability mass function, $P(D_{\sf HARQ}=d$)';
        box on;     set(gca,'FontSize',20) 
        xlhand = get(gca,'xlabel'); xlabel(xlab,'Interpreter','latex'); set(xlhand,'fontsize',20) 
        ylhand = get(gca,'ylabel'); ylabel(ylab,'Interpreter','latex'); set(ylhand,'fontsize',20)
        xaxis = T;    yaxis = 0.1;

        if K==1
            text(xaxis,yaxis,['$\sum\limits_{n=0}^{5T} P(D_{\sf HARQ}=n)$=' num2str( sum(PMF_delay(PMF_delay>0)) ) ', memoryless' ', $\epsilon$=' num2str(epsf) ', $T$=' num2str(T) ', $k$=' num2str(k)],'fontsize',20,'Interpreter','latex');
        else
            text(xaxis,yaxis,['$\sum\limits_{n=0}^{5T} P(D_{\sf HARQ}=n)$=' num2str( sum(PMF_delay(PMF_delay>0)) ) ', Gilbert-Elliott' ', $\epsilon$=' num2str(epsf) ', $T$=' num2str(T) ', $k$=' num2str(k)],'fontsize',20,'Interpreter','latex');
        end
    else %Ztransform==0                 
        %Numerical integration (assume a continuous delay distribution)
        %z = exp(1i*t)  
        f = @(t,x) 1/(2*pi) * exp(-1i*t.*x).* phiD(exp(1i*t));
        %f_phiD = integral(@(t) f(t,x),-Inf,Inf);
        tmin = -pi;
        tmax = pi;
        tno = 10^5;
        t_range = linspace(tmin,tmax,tno);
        t_delta = t_range(2)- t_range(1);
        n_min = k;
        n_max = 10*T;
        no = n_max-n_min+1;
        lx = 1000*no;
        x_range = linspace(n_min,n_max,lx);
        x_delta = x_range(2)-x_range(1);
        f_phiD = zeros(1,lx);    
        for t = t_range
            f_phiD = f_phiD + f(t,x_range)*t_delta;
        end  
        PMF_delay = zeros(1,no);
        for count = 1:no
            n = n_min+count-1
            indices = find(x_range>n-1 & x_range<=n);
            PMF_delay(count) = real(sum(f_phiD(indices))*x_delta);
        end
        figure    
        plot(n_min:n_max,PMF_delay,'k','linewidth',2) %plot(x_range,real(f_phiD),'k','linewidth',2)
        xlab = 'Value of delay, d'; 
        ylab = 'Probability mass function, $P(D_{\sf HARQ}=d$)';
        box on;     set(gca,'FontSize',20) 
        xlhand = get(gca,'xlabel'); xlabel(xlab,'Interpreter','latex'); set(xlhand,'fontsize',20) 
        ylhand = get(gca,'ylabel'); ylabel(ylab,'Interpreter','latex'); set(ylhand,'fontsize',20)
        xaxis = T;    yaxis = 0.1;

        if K==1
            text(xaxis,yaxis,['$\sum\limits_{n=0}^{10T} P(D_{\sf HARQ}=n)$=' num2str( sum(PMF_delay(PMF_delay>0)) ) ', memoryless' ', $\epsilon$=' num2str(epsf) ', $T$=' num2str(T) ', $k$=' num2str(k)],'fontsize',20,'Interpreter','latex');
        else
            text(xaxis,yaxis,['$\sum\limits_{n=0}^{10T} P(D_{\sf HARQ}=n)$=' num2str( sum(PMF_delay(PMF_delay>0)) ) ', Gilbert-Elliott' ', $\epsilon$=' num2str(epsf) ', $T$=' num2str(T) ', $k$=' num2str(k)],'fontsize',20,'Interpreter','latex');
        end
    end
    %CCDF:P(D>d)=1-P(D<=d)
    CCDF_delay = zeros(1,n_min+no-1);
    CCDF_delay(1:n_min-1) = 1;
    for count = 1:no
        n = n_min+count-1;
        CCDF_delay(n) = 1-sum(PMF_delay(1:count));
    end
    figure    
    plot(1:n_max,CCDF_delay,'k','linewidth',2) 
    xlab = 'Value of delay, d'; 
    ylab = 'CCDF, $P(D_{\sf HARQ}>d$)';
    box on;     set(gca,'FontSize',20) 
    xlhand = get(gca,'xlabel'); xlabel(xlab,'Interpreter','latex'); set(xlhand,'fontsize',20) 
    ylhand = get(gca,'ylabel'); ylabel(ylab,'Interpreter','latex'); set(ylhand,'fontsize',20)
    xaxis = T;    yaxis = 0.1;
    if K==1
        text(xaxis,yaxis,['$\sum\limits_{n=0}^{5T} P(D_{\sf HARQ}=n)$=' num2str( sum(PMF_delay(PMF_delay>0)) ) ', memoryless' ', $\epsilon$=' num2str(epsf) ', $T$=' num2str(T) ', $k$=' num2str(k)],'fontsize',20,'Interpreter','latex');
    else
        text(xaxis,yaxis,['$\sum\limits_{n=0}^{5T} P(D_{\sf HARQ}=n)$=' num2str( sum(PMF_delay(PMF_delay>0)) ) ', Gilbert-Elliott' ', $\epsilon$=' num2str(epsf) ', $T$=' num2str(T) ', $k$=' num2str(k)],'fontsize',20,'Interpreter','latex');
    end    

end

