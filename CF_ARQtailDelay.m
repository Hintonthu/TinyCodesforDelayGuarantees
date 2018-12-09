function [phiD, PMF_delay, CCDF_delay] = CF_ARQtailDelay(k,K,T,P_kron,pi_I_kron,P00,P01,P10,P11,P0x,Px0,P1x,Px1,epsf,M,Ztransform)

    syms z

    RTT = k+M-1;       %k = RTT-M+1; %time needed to get the first feedback after the last packet is transmitted
    d = T-RTT;         %time to timeout    
   
    t = T; 
    if K~=1
        I = eye(2*K,2*K);
        onevec = ones(2*K,1);
        zeromat = zeros(2*K,2*K);
    else %K==1
        I = 1;
        onevec = 1;
        zeromat = 0;
    end
    I2 = I; 
    
    
    %b*inv(A) = b/A         %inv(A)*b = A\b
    
    P = P_kron;
    
    % ANALYTICAL EXPRESSION FOR DELAY   
    P10C = P10*P10+P10*P01+P01*P10;
    P11C = P11*P11+P11*P10+P10*P11;
    
    P1xC = P10C+P11C;

    
    P00C = P00*P00;
    P01C = P01*P01+P01*P00+P00*P01;
    
    P0xC = P00C+P01C;

    
    P0 = P00*P10+P10*P00...
       + P11*P01+P01*P11... + P01*P00+P00*P01...
       + P00*P11+P11*P00;
    %P10C+P11C+P00C+P01C+P0 = I2
    %equivalently P1XC+P0XC+P0 = I2
    %equivalently PX0C+PX1C+P0 = I2
    
    Px0C = P00C+P10C; %Px0*Px0;
    Px1C = P01C+P11C; %Px1*Px1;
    
    
    
    P1x_D_1 = @(z) z^(T-d)*P^(T-d)*(z*P10C+z*P11C*z^d*P^d);
    P1x_D_2 = @(z) z^(T-d-1)*P^(T-d-1)*(z*P10+z*P11*z^(d+1)*P^(d+1));
    %STATE C2 & self loop: Px0C -> (I2-Px1C) because in my definition, Px0C+Px1C \neq I2
    PhiD = @(z) z^k*P^k/(I2-P1x_D_1(z))*(z*P00C+z*P01C/(I2-z*Px1C)*z*(I2-Px1C)+z*P0/(I2-P1x_D_2(z))*(z*P00+z*P01/(I2-z*Px1)*z*Px0));
    phiD = @(z) (pi_I_kron * PhiD(z) * onevec)/(pi_I_kron * onevec);
    
    if Ztransform==1
        %Z-transform of x[n]: sum_n=-inf^inf x[n]z^-n
        %If x[n]=P(X=n), and z=z^-1, this gives the probability-generating function (PGF) for X: E[z^X].
        n_min = 5*T+1;%k+M-1;
        n_max = 10*T;%5*T;
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
        ylab = 'Probability mass function, $P(D_{\sf CF-ARQ}=d$)';
        box on;     set(gca,'FontSize',20) 
        xlhand = get(gca,'xlabel'); xlabel(xlab,'Interpreter','latex'); set(xlhand,'fontsize',20) 
        ylhand = get(gca,'ylabel'); ylabel(ylab,'Interpreter','latex'); set(ylhand,'fontsize',20)
        xaxis = T;    yaxis = 0.1;

        if K==1
            text(xaxis,yaxis,['$\sum\limits_{n=0}^{5T} P(D_{\sf CF-ARQ}=n)$=' num2str( sum(PMF_delay(PMF_delay>0)) ) ', memoryless' ', $\epsilon$=' num2str(epsf) ', $T$=' num2str(T) ', $k$=' num2str(k)],'fontsize',20,'Interpreter','latex');
        else
            text(xaxis,yaxis,['$\sum\limits_{n=0}^{5T} P(D_{\sf CF-ARQ}=n)$=' num2str( sum(PMF_delay(PMF_delay>0)) ) ', Gilbert-Elliott' ', $\epsilon$=' num2str(epsf) ', $T$=' num2str(T) ', $k$=' num2str(k)],'fontsize',20,'Interpreter','latex');
        end
        
    else %Ztransform==0
        %Numerical integration (assume a continuous delay distribution)     
        f = @(t,x) 1/(2*pi) * exp(-1i*t.*x).* phiD(exp(1i*t));
        %f_phiD = integral(@(t) f(t,x),-Inf,Inf);
        tmin = -pi;
        tmax = pi;
        tno = 10^5;
        t_range = linspace(tmin,tmax,tno);
        t_delta = t_range(2)- t_range(1);
        n_min = 1;%k+M-1;
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
            indices = x_range>n-1 & x_range<=n;
            PMF_delay(count) = real(sum(f_phiD(indices)))*x_delta;
        end
        figure    
        plot(n_min:n_max,PMF_delay,'k','linewidth',2)
        xlab = 'Value of delay, d'; 
        ylab = 'Probability mass function, $P(D_{\sf CF-ARQ}=d$)';
        box on;     set(gca,'FontSize',20) 
        xlhand = get(gca,'xlabel'); xlabel(xlab,'Interpreter','latex'); set(xlhand,'fontsize',20) 
        ylhand = get(gca,'ylabel'); ylabel(ylab,'Interpreter','latex'); set(ylhand,'fontsize',20)
        xaxis = T;    yaxis = 0.1;

        if K==1
            text(xaxis,yaxis,['$\sum\limits_{n=0}^{10T} P(D_{\sf CF-ARQ}=n)$=' num2str( sum(PMF_delay(k+M-1:count)) ) ', memoryless' ', $\epsilon$=' num2str(epsf) ', $T$=' num2str(T) ', $k$=' num2str(k)],'fontsize',20,'Interpreter','latex');
        else
            text(xaxis,yaxis,['$\sum\limits_{n=0}^{10T} P(D_{\sf CF-ARQ}=n)$=' num2str( sum(PMF_delay(k+M-1:count)) ) ', Gilbert-Elliott' ', $\epsilon$=' num2str(epsf) ', $T$=' num2str(T) ', $k$=' num2str(k)],'fontsize',20,'Interpreter','latex');
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
    ylab = 'CCDF, $P(D_{\sf CF-ARQ}>d$)';
    box on;     set(gca,'FontSize',20) 
    xlhand = get(gca,'xlabel'); xlabel(xlab,'Interpreter','latex'); set(xlhand,'fontsize',20) 
    ylhand = get(gca,'ylabel'); ylabel(ylab,'Interpreter','latex'); set(ylhand,'fontsize',20)
    xaxis = T;    yaxis = 0.1;
    if K==1
        text(xaxis,yaxis,['$\sum\limits_{n=0}^{5T} P(D_{\sf CF-ARQ}=n)$=' num2str( sum(PMF_delay(PMF_delay>0)) ) ', memoryless' ', $\epsilon$=' num2str(epsf) ', $T$=' num2str(T) ', $k$=' num2str(k)],'fontsize',20,'Interpreter','latex');
    else
        text(xaxis,yaxis,['$\sum\limits_{n=0}^{5T} P(D_{\sf CF-ARQ}=n)$=' num2str( sum(PMF_delay(PMF_delay>0)) ) ', Gilbert-Elliott' ', $\epsilon$=' num2str(epsf) ', $T$=' num2str(T) ', $k$=' num2str(k)],'fontsize',20,'Interpreter','latex');
    end
end

