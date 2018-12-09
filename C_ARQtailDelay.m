function [phiD, PMF_delay, CCDF_delay] = C_ARQtailDelay(k,K,T,P_kron,pi_I_kron,P00,P01,P10,P11,P0x,Px0,P1x,Px1,epsf,M,Ztransform);
    syms z
    d = T-(k+M-1);          %time to timeout    
     
    if K~=1
        I = eye(2*K,2*K);
        onevec = ones(2*K,1);
        zeromat = zeros(2*K,2*K);
    else %K==1
        I = 1;
        onevec = 1;
        zeromat = 0;
    end
    I2 = I; %[I, zeromat;zeromat, I];
    
    
    %b*inv(A) = b/A         %inv(A)*b = A\b
    
    P = P_kron;
    
    % ANALYTICAL EXPRESSION FOR DELAY     
    P0x = P00+P01;
    P1x = P11+P10;
    P000 = P0x*P00;
    P001 = P0x*P01;
    P010 = P0x*P10;
    P011 = P0x*P11;
    P100 = P1x*P00;
    P101 = P1x*P01;
    P110 = P1x*P10;
    P111 = P1x*P11;
    Pxx0 = P000+P010+P100+P110;
    Pxx1 = P001+P011+P101+P111;
    
    % Functions for delay analysis
    f1 = @(z) (z*P10+z*P11*z^(T-k)*P^(T-k))*z^(k-1)*P^(k-1);
    f2 = @(z) (z*P110+z*P111*z^(T-k)*P^(T-k))*z^k*P^k;
    f3 = @(z) z*P11*z^T*P^T;
    
    A1 = @(z) z*(P100+P010)+z*(P011+P101)*(z^T*P^T+I2/(I2-f3(z))-I2)*z*P10;              
    A2 = @(z) z*P00+z*P01/(I2-z*Px1)*z*Px0;
    A3 = @(z) z*P000+z*P001/(I2-z*Pxx1)*z*Pxx0;
    A4 = @(z) z*(P011+P101)*(z^T*P^T+I2/(I2-f3(z))-I2)*(z*P00+z*P01/(I2-z*Px1)*z*Px0);
    
    PhiD = @(z) z^k*P^k/(I2-f2(z))*( A1(z)*I2/(I2-f1(z))*A2(z)+A3(z)+A4(z));
    
    %or equivalently
    %B00 = @(z) z*P00+z*P01/(I2-z*Px1)*z*Px0;
    %B000 = @(z) z*P000+z*P001/(I2-z*Pxx1)*z*Pxx0;
    %PhiD = @(z) z^k*P^k/(I2-f2(z))*[[z*(P100+P010)+z*(P011+P101)*(z^T*P^T+I2/(I2-f3(z))-I2)*z*P10]*I2/(I2-f1(z))*B00(z)+B000(z)+z*(P011+P101)*(z^T*P^T+I2/(I2-f3(z))-I2)*B00(z)];
      
    phiD = @(z) (pi_I_kron * PhiD(z) * onevec)/(pi_I_kron * onevec);
    
    if Ztransform==1
        %Z-transform of x[n]: sum_n=-inf^inf x[n]z^-n
        %If x[n]=P(X=n), and z=z^-1, this gives the probability-generating function (PGF) for X: E[z^X].
        n_min = k+M-1;
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
        ylab = 'Probability mass function, $P(D_{\sf C-ARQ}=d$)';
        box on;     set(gca,'FontSize',20) 
        xlhand = get(gca,'xlabel'); xlabel(xlab,'Interpreter','latex'); set(xlhand,'fontsize',20) 
        ylhand = get(gca,'ylabel'); ylabel(ylab,'Interpreter','latex'); set(ylhand,'fontsize',20)
        xaxis = T;    yaxis = 0.1;

        if K==1
            text(xaxis,yaxis,['$\sum\limits_{n=0}^{5T} P(D_{\sf C-ARQ}=n)$=' num2str( sum(PMF_delay(PMF_delay>0)) ) ', memoryless' ', $\epsilon$=' num2str(epsf) ', $T$=' num2str(T) ', $k$=' num2str(k)],'fontsize',20,'Interpreter','latex');
        else
            text(xaxis,yaxis,['$\sum\limits_{n=0}^{5T} P(D_{\sf C-ARQ}=n)$=' num2str( sum(PMF_delay(PMF_delay>0)) ) ', Gilbert-Elliott' ', $\epsilon$=' num2str(epsf) ', $T$=' num2str(T) ', $k$=' num2str(k)],'fontsize',20,'Interpreter','latex');
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
            indices = x_range>n-2 & x_range<=n-1;
            PMF_delay(count) = real(sum(f_phiD(indices)))*x_delta;
        end
        figure    
        plot(n_min:n_max,PMF_delay,'k','linewidth',2)
        xlab = 'Value of delay, d'; 
        ylab = 'Probability mass function, $P(D_{\sf C-ARQ}=d$)';
        box on;     set(gca,'FontSize',20) 
        xlhand = get(gca,'xlabel'); xlabel(xlab,'Interpreter','latex'); set(xlhand,'fontsize',20) 
        ylhand = get(gca,'ylabel'); ylabel(ylab,'Interpreter','latex'); set(ylhand,'fontsize',20)
        xaxis = T;    yaxis = 0.1;

        if K==1
            text(xaxis,yaxis,['$\sum\limits_{n=0}^{10T} P(D_{\sf C-ARQ}=n)$=' num2str( sum(PMF_delay(k+M-1:count)) ) ', memoryless' ', $\epsilon$=' num2str(epsf) ', $T$=' num2str(T) ', $k$=' num2str(k)],'fontsize',20,'Interpreter','latex');
        else
            text(xaxis,yaxis,['$\sum\limits_{n=0}^{10T} P(D_{\sf C-ARQ}=n)$=' num2str( sum(PMF_delay(k+M-1:count)) ) ', Gilbert-Elliott' ', $\epsilon$=' num2str(epsf) ', $T$=' num2str(T) ', $k$=' num2str(k)],'fontsize',20,'Interpreter','latex');
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
    ylab = 'CCDF, $P(D_{\sf C-ARQ}>d$)';
    box on;     set(gca,'FontSize',20) 
    xlhand = get(gca,'xlabel'); xlabel(xlab,'Interpreter','latex'); set(xlhand,'fontsize',20) 
    ylhand = get(gca,'ylabel'); ylabel(ylab,'Interpreter','latex'); set(ylhand,'fontsize',20)
    xaxis = T;    yaxis = 0.1;
    if K==1
        text(xaxis,yaxis,['$\sum\limits_{n=0}^{5T} P(D_{\sf C-ARQ}=n)$=' num2str( sum(PMF_delay(PMF_delay>0)) ) ', memoryless' ', $\epsilon$=' num2str(epsf) ', $T$=' num2str(T) ', $k$=' num2str(k)],'fontsize',20,'Interpreter','latex');
    else
        text(xaxis,yaxis,['$\sum\limits_{n=0}^{5T} P(D_{\sf C-ARQ}=n)$=' num2str( sum(PMF_delay(PMF_delay>0)) ) ', Gilbert-Elliott' ', $\epsilon$=' num2str(epsf) ', $T$=' num2str(T) ', $k$=' num2str(k)],'fontsize',20,'Interpreter','latex');
    end       
end

