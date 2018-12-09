function [phiD, PMF_delay, CCDF_delay] = UncodedARQtailDelay(k,K,T,P_kron,pi_I_kron,P00,P01,P10,P11,P0x,Px0,P1x,Px1,epsf,Ztransform)
    
    syms z

    d = T-k;               %time to timeout  
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

    A = @(z) I-z^k*P10*P_kron^(k-1)-z^T*P11*P_kron^(T-1);
    %B = @(z) I-z*P10*P_kron^(k-1)-z*P11*P_kron^(t-1);
    %Btau = @(z) I-z*P10*P_kron^(k-1)-z*P11*P_kron^(T-1);
    %C = I-Px1;

    PhiD = @(z) z^(k-1)*P_kron^(k-1)/A(z)*(z*P00+z^2*(P01/(I-z*Px1))*Px0);    

    phiD = @(z) (pi_I_kron * PhiD(z) * onevec)/(pi_I_kron * onevec);
    
    
    if Ztransform == 1
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
        ylab = 'Probability mass function, $P(D_{\sf ARQ}=d$)';
        box on;     set(gca,'FontSize',20) 
        xlhand = get(gca,'xlabel'); xlabel(xlab,'Interpreter','latex'); set(xlhand,'fontsize',20) 
        ylhand = get(gca,'ylabel'); ylabel(ylab,'Interpreter','latex'); set(ylhand,'fontsize',20)
        xaxis = T;    yaxis = 0.1;
        text(xaxis,yaxis,['$\sum\limits_{n=0}^{5T} P(D_{\sf ARQ}=n)$=' num2str( sum(PMF_delay(PMF_delay>0)) ) ', memoryless' ', $\epsilon$=' num2str(epsf) ', $T$=' num2str(T) ', $k$=' num2str(k)],'fontsize',20,'Interpreter','latex');
        if K==1
            text(xaxis,yaxis,['$\sum\limits_{n=0}^{5T} P(D_{\sf ARQ}=n)$=' num2str( sum(PMF_delay(PMF_delay>0)) ) ', memoryless' ', $\epsilon$=' num2str(epsf) ', $T$=' num2str(T) ', $k$=' num2str(k)],'fontsize',20,'Interpreter','latex');
        else
            text(xaxis,yaxis,['$\sum\limits_{n=0}^{5T} P(D_{\sf ARQ}=n)$=' num2str( sum(PMF_delay(PMF_delay>0)) ) ', Gilbert-Elliott' ', $\epsilon$=' num2str(epsf) ', $T$=' num2str(T) ', $k$=' num2str(k)],'fontsize',20,'Interpreter','latex');
        end
        
    else %Ztransform == 0
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
        n_max = 5*T;
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
        ylab = 'Probability mass function, $P(D_{\sf ARQ}=d$)';
        box on;     set(gca,'FontSize',20) 
        xlhand = get(gca,'xlabel'); xlabel(xlab,'Interpreter','latex'); set(xlhand,'fontsize',20) 
        ylhand = get(gca,'ylabel'); ylabel(ylab,'Interpreter','latex'); set(ylhand,'fontsize',20)
        xaxis = T;    yaxis = 0.1;
        text(xaxis,yaxis,['$\sum\limits_{n=0}^{10T} P(D_{\sf ARQ}=n)$=' num2str( sum(PMF_delay(PMF_delay>0)) ) ', Gilbert-Elliott' ', $\epsilon$=' num2str(epsf) ', $T$=' num2str(T) ', $k$=' num2str(k)],'fontsize',20,'Interpreter','latex');
        if K==1
            text(xaxis,yaxis,['$\sum\limits_{n=0}^{5T} P(D_{\sf ARQ}=n)$=' num2str( sum(PMF_delay(PMF_delay>0)) ) ', memoryless' ', $\epsilon$=' num2str(epsf) ', $T$=' num2str(T) ', $k$=' num2str(k)],'fontsize',20,'Interpreter','latex');
        else
            text(xaxis,yaxis,['$\sum\limits_{n=0}^{5T} P(D_{\sf ARQ}=n)$=' num2str( sum(PMF_delay(PMF_delay>0)) ) ', Gilbert-Elliott' ', $\epsilon$=' num2str(epsf) ', $T$=' num2str(T) ', $k$=' num2str(k)],'fontsize',20,'Interpreter','latex');
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
    ylab = 'CCDF, $P(D_{\sf ARQ}>d$)';
    box on;     set(gca,'FontSize',20) 
    xlhand = get(gca,'xlabel'); xlabel(xlab,'Interpreter','latex'); set(xlhand,'fontsize',20) 
    ylhand = get(gca,'ylabel'); ylabel(ylab,'Interpreter','latex'); set(ylhand,'fontsize',20)
    xaxis = T;    yaxis = 0.1;
    if K==1
        text(xaxis,yaxis,['$\sum\limits_{n=0}^{5T} P(D_{\sf ARQ}=n)$=' num2str( sum(PMF_delay(PMF_delay>0)) ) ', memoryless' ', $\epsilon$=' num2str(epsf) ', $T$=' num2str(T) ', $k$=' num2str(k)],'fontsize',20,'Interpreter','latex');
    else
        text(xaxis,yaxis,['$\sum\limits_{n=0}^{5T} P(D_{\sf ARQ}=n)$=' num2str( sum(PMF_delay(PMF_delay>0)) ) ', Gilbert-Elliott' ', $\epsilon$=' num2str(epsf) ', $T$=' num2str(T) ', $k$=' num2str(k)],'fontsize',20,'Interpreter','latex');
    end    

end

