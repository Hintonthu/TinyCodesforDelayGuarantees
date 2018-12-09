function [PhiD, phiD, meanDelay, varDelay, Throughput] = CodedPhiD(z,k,K,T,P_kron,pi_I_kron,P00,P01,P10,P11,Px0,Px1,M,N,symbolic)
%Matrix-generating function of delay, coding, new model

    %M = 2 packets transmitted
    %N = M = 2 DoFs required
    if symbolic %for memoryless model
      
        syms e T k z
        P00 = (1-e)^2;
        P01 = (1-e)*e;
        P10 = e*(1-e);
        P11 = e^2;

        P0x = P00+P01;
        P1x = P10+P11;
        Px0 = P00+P10;
        Px1 = P01+P11; 
        
        e_vec = e;
        q = 0;
        
        P = 1;
        P_kron = 1;        
        K = 1;
        
        P1 = P * diag(e_vec);       %error probability matrix of HMM
        P0 = P * diag(1-e_vec);     %success probability matrix of HMM

        pivec = 1;
             
        pi_kron = kron(pivec,pivec);
        pi_I = pivec*P0; 
        pi_I_kron = kron(pi_I,pi_I); 
    end
    
    
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
    f1 = (z*P10+z*P11*z^(T-k)*P^(T-k))*z^(k-1)*P^(k-1);
    f2 = (z*P110+z*P111*z^(T-k)*P^(T-k))*z^k*P^k;
    f3 = z*P11*z^T*P^T;
    
    A1 = z*(P100+P010)+z*(P011+P101)*(z^T*P^T+I2/(I2-f3))*z*P10;              
    A2 = z*P00+z*P01/(I2-z*Px1)*z*Px0;
    A3 = z*P000+z*P001/(I2-z*Pxx1)*z*Pxx0;
    A4 = z*(P011+P101)*(z^T*P^T+I2/(I2-f3))*(z*P00+z*P01/(I2-z*Px1)*z*Px0);
    
    PhiD = z^k*P^k/(I2-f2)*( A1*I2/(I2-f1)*A2+A3+A4);
      
    phiD = (pi_I_kron * PhiD * onevec)/(pi_I_kron * onevec);
    
    %Derivative of phiD: evaluate at z=1 to obtain mean delay
    
    derivativef1 = k*P10*z^(k-1)*P^k+T*P11*z^(T-1)*P^T;
    derivativef2 = (k+1)*P110*z^k*P^k+(T+1)*P111*z^T*P^T;
    derivativef3 = (T+1)*P11*z^T*P^T;
    
    derivativeA1 = (P100+P010)+(P011+P101)*((T+2)*z^(T+1)*P^T+2*z*I2/(I2-f3)+z^2*I2/(I2-f3)*derivativef3/(I2-f3))*P10;
    derivativeA2 = P00+2*z*P01/(I2-z*Px1)*Px0+z^2*P01/(I2-z*Px1)*Px1/(I2-z*Px1)*Px0;
    derivativeA3 = P000+2*z*P001/(I2-z*Pxx1)*Pxx0+z^2*P001/(I2-z*Pxx1)*Pxx1/(I2-z*Pxx1)*Pxx0;
    derivativeA4 = (P011+P101)*(z^T*P^T+I2/(I2-f3))*(z*P00+z*P01/(I2-z*Px1)*z*Px0)...
                 + z*(P011+P101)*(T*z^(T-1)*P^T+I2/(I2-f3)*derivativef3/(I2-f3))*(z*P00+z*P01/(I2-z*Px1)*z*Px0)...
                 + z*(P011+P101)*(z^T*P^T+I2/(I2-f3))*(P00+2*z*P01/(I2-z*Px1)*Px0+z^2*P01/(I2-z*Px1)*Px1/(I2-z*Px1)*Px0);
                    
    derivativePhiD = k*z^(k-1)*P^k/(I2-f2)*( A1*I2/(I2-f1)*A2+A3+A4)...
                   + z^k*P^k/(I2-f2)*derivativef2/(I2-f2)*( A1*I2/(I2-f1)*A2+A3+A4)...
                   + z^k*P^k/(I2-f2)*( derivativeA1*I2/(I2-f1)*A2+A1*I2/(I2-f1)*derivativef1/(I2-f1)*A2+A1*I2/(I2-f1)*derivativeA2+derivativeA3+derivativeA4);
    

    %Mean delay per M packets (bucket)           
    meanDelay = (pi_I_kron * derivativePhiD * onevec)/(pi_I_kron * onevec); 
    
    
    % ANALYTICAL EXPRESSION FOR THROUGHPUT 
    
    % Functions for throughput analysis
    g1 = (P10+P11*P^(T-k))*z*P^(k-1);
    g2 = (P110+P111*P^(T-k))*z*P^k;
    g3 = z*P11*P^T;
    
     if symbolic
        sum1 = (1-Px1^T)/(1-Px1);
        sum2 = (1-Px1^(d+1))/(1-Px1);
        sum3 = (1-Pxx1^(d+1))/(1-Pxx1);
        sum4 = (1-Pxx1^T)/(1-Pxx1);
     else
        sum1 = I2;
        sum2 = I2;
        sum3 = I2;
        sum4 = I2;
        for i = 1:T-1
            sum1 = sum1 + Px1^i;           
            if i<=d
                sum2 = sum2 + Px1^i;   %1st layer where 2DoFs required   
                sum3 = sum3 + Pxx1^i;   
            end
            sum4 = sum4 + Pxx1^i;   
        end
     end
    
    B1 = ((P100+P010)+(P011+P101)*(z*P^T+I2/(I2-g3))*P10)/(I2-g1)*((P00+P01*sum2*Px0)+(P01*Px1^(T-k)/(I2-z*Px1^T)*z*sum1*Px0));
    B2 = P000+P001*sum3*Pxx0;
    B3 = P001*Pxx1^(T-k)/(I2-z*Pxx1^T)*z*sum4*Pxx0;
    B4 = (P011+P101)*(z*P^T+I2/(I2-g3))*((P00+P01*sum2*Px0)+(P01*Px1^(T-k)/(I2-z*Px1^T)*z*sum1*Px0));
 
    
    Phitau = z^2*P^k/(I2-g2)*(B1+B2+B3+B4);    
    
    derivativeg1 = (P10+P11*P^(T-k))*P^(k-1);
    derivativeg2 = (P110+P111*P^(T-k))*P^k;
    derivativeg3 = P11*P^T;
    
    derivativeB1 = ((P011+P101)*(P^T+I2/(I2-g3)*derivativeg3/(I2-g3))*P10)/(I2-g1)*((P00+P01*sum2*Px0)+(P01*Px1^(T-k)/(I2-z*Px1^T)*z*sum1*Px0))...
                 + ((P100+P010)+(P011+P101)*(z*P^T+I2/(I2-g3))*P10)/(I2-g1)*derivativeg1/(I2-g1)*((P00+P01*sum2*Px0)+(P01*Px1^(T-k)/(I2-z*Px1^T)*z*sum1*Px0))...
                 + ((P100+P010)+(P011+P101)*(z*P^T+I2/(I2-g3))*P10)/(I2-g1)*(P01*Px1^(T-k)/(I2-z*Px1^T)*Px1^T/(I2-z*Px1^T)*z*sum1*Px0+P01*Px1^(T-k)/(I2-z*Px1^T)*sum1*Px0);
    derivativeB2 = zeromat;
    derivativeB3 = P001*Pxx1^(T-k)/(I2-z*Pxx1^T)*Pxx1^T/(I2-z*Pxx1^T)*z*sum4*Pxx0+P001*Pxx1^(T-k)/(I2-z*Pxx1^T)*sum4*Pxx0;
    derivativeB4 = (P011+P101)*(P^T+I2/(I2-g3)*derivativeg3/(I2-g3))*((P00+P01*sum2*Px0)+(P01*Px1^(T-k)/(I2-z*Px1^T)*z*sum1*Px0))...
                 + (P011+P101)*(z*P^T+I2/(I2-g3))*(P01*Px1^(T-k)/(I2-z*Px1^T)*Px1^T/(I2-z*Px1^T)*z*sum1*Px0+P01*Px1^(T-k)/(I2-z*Px1^T)*sum1*Px0);
        
    
    derivativePhitau = 2*z*P^k/(I2-g2)*(B1+B2+B3+B4)...
                     + z^2*P^k/(I2-g2)*derivativeg2/(I2-g2)*(B1+B2+B3+B4)...
                     + z^2*P^k/(I2-g2)*(derivativeB1+derivativeB2+derivativeB3+derivativeB4);
                                 
     
    %Total transmission time for M packets                             
    transmissionTime = (pi_I_kron * derivativePhitau * onevec)/(pi_I_kron * onevec); 
    
    Throughput = 1/transmissionTime;    
    
    
    %VARIABILITY OF DELAY: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    secondderivativef1 = k*(k-1)*P10*z^(k-2)*P^k+T*(T-1)*P11*z^(T-2)*P^T;
    secondderivativef2 = (k+1)*k*P110*z^(k-1)*P^k+(T+1)*T*P111*z^(T-1)*P^T;
    secondderivativef3 = (T+1)*T*P11*z^(T-1)*P^T;
    
    secondderivativeA1 = (P011+P101)*((T+2)*(T+1)*z^T*P^T+2*I2/(I2-f3)+2*z*I2/(I2-f3)*derivativef3/(I2-f3)+2*z*I2/(I2-f3)*derivativef3/(I2-f3)+z^2*I2/(I2-f3)*derivativef3/(I2-f3)*derivativef3/(I2-f3)+z^2*I2/(I2-f3)*secondderivativef3/(I2-f3)+z^2*I2/(I2-f3)*derivativef3/(I2-f3)*derivativef3/(I2-f3))*P10;
    secondderivativeA2 = 2*P01/(I2-z*Px1)*Px0+2*z*P01/(I2-z*Px1)*Px1/(I2-z*Px1)*Px0+z^2*P01/(I2-z*Px1)*Px1/(I2-z*Px1)*Px1/(I2-z*Px1)*Px0 +z^2*P01/(I2-z*Px1)*Px1/(I2-z*Px1)*Px1/(I2-z*Px1)*Px0; 
    secondderivativeA3 = 2*P001/(I2-z*Pxx1)*Pxx0+ 2*z*P001/(I2-z*Pxx1)*Pxx1/(I2-z*Pxx1)*Pxx0 +2*z*P001/(I2-z*Pxx1)*Pxx1/(I2-z*Pxx1)*Pxx0 +z^2*P001/(I2-z*Pxx1)*Pxx1 /(I2-z*Pxx1)*Pxx1/(I2-z*Pxx1)*Pxx0 +z^2*P001/(I2-z*Pxx1)*Pxx1/(I2-z*Pxx1)*Pxx1/(I2-z*Pxx1)*Pxx0;
    secondderivativeA4 = (P011+P101)*(T*z^(T-1)*P^T+I2/(I2-f3)*derivativef3/(I2-f3))*(z*P00+z*P01/(I2-z*Px1)*z*Px0)+(P011+P101)*(z^T*P^T+I2/(I2-f3))*(P00+2*z*P01/(I2-z*Px1)*Px0+z*P01/(I2-z*Px1)*Px1/(I2-z*Px1)*z*Px0)...
                 + (P011+P101)*(T*z^(T-1)*P^T+I2/(I2-f3)*derivativef3/(I2-f3))*(z*P00+z*P01/(I2-z*Px1)*z*Px0)+ z*(P011+P101)*(T*(T-1)*z^(T-2)*P^T+I2/(I2-f3)*derivativef3/(I2-f3)*derivativef3/(I2-f3)+I2/(I2-f3)*secondderivativef3/(I2-f3)+I2/(I2-f3)*derivativef3/(I2-f3)*derivativef3/(I2-f3))*(z*P00+z*P01/(I2-z*Px1)*z*Px0)+z*(P011+P101)*(T*z^(T-1)*P^T+I2/(I2-f3)*derivativef3/(I2-f3))*(P00+2*z*P01/(I2-z*Px1)*Px0+z*P01/(I2-z*Px1)*Px1/(I2-z*Px1)*z*Px0)...
                 + (P011+P101)*(z^T*P^T+I2/(I2-f3))*(P00+2*z*P01/(I2-z*Px1)*Px0+z^2*P01/(I2-z*Px1)*Px1/(I2-z*Px1)*Px0) + z*(P011+P101)*(T*z^(T-1)*P^T+I2/(I2-f3)*derivativef3/(I2-f3))*(P00+2*z*P01/(I2-z*Px1)*Px0+z^2*P01/(I2-z*Px1)*Px1/(I2-z*Px1)*Px0) + z*(P011+P101)*(z^T*P^T+I2/(I2-f3))*(2*P01/(I2-z*Px1)*Px0+2*z*P01/(I2-z*Px1)*Px1/(I2-z*Px1)*Px0+2*z*P01/(I2-z*Px1)*Px1/(I2-z*Px1)*Px0 +z^2*P01/(I2-z*Px1)*Px1/(I2-z*Px1)*Px1/(I2-z*Px1)*Px0 +z^2*P01/(I2-z*Px1)*Px1/(I2-z*Px1)*Px1/(I2-z*Px1)*Px0);
    
     
    %b*inv(A) = b/A         %inv(A)*b = A\b
    
    %Second derivative of phiD: evaluate at z=1 to obtain the second moment of delay
    secondderivativePhiD = k*(k-1)*z^(k-2)*P^k/(I2-f2)*( A1*I2/(I2-f1)*A2+A3+A4)+k*z^(k-1)*P^k/(I2-f2)*derivativef2/(I2-f2)*( A1*I2/(I2-f1)*A2+A3+A4)+k*z^(k-1)*P^k/(I2-f2)*( derivativeA1*I2/(I2-f1)*A2+A1*I2/(I2-f1)*derivativef1/(I2-f1)*A2+A1*I2/(I2-f1)*derivativeA2+derivativeA3+derivativeA4)...
                   + k*z^(k-1)*P^k/(I2-f2)*derivativef2/(I2-f2)*( A1*I2/(I2-f1)*A2+A3+A4)+z^k*P^k/(I2-f2)*derivativef2/(I2-f2)*derivativef2/(I2-f2)*( A1*I2/(I2-f1)*A2+A3+A4)+z^k*P^k/(I2-f2)*derivativef2/(I2-f2)*( A1*I2/(I2-f1)*A2+A3+A4)+z^k*P^k/(I2-f2)*secondderivativef2/(I2-f2)*( A1*I2/(I2-f1)*A2+A3+A4)+z^k*P^k/(I2-f2)*derivativef2/(I2-f2)*derivativef2/(I2-f2)*( A1*I2/(I2-f1)*A2+A3+A4)+z^k*P^k/(I2-f2)*derivativef2/(I2-f2)*( derivativeA1*I2/(I2-f1)*A2+A1*I2/(I2-f1)*derivativef1/(I2-f1)*A2+A1*I2/(I2-f1)*derivativeA2+derivativeA3+derivativeA4)...
                   + k*z^(k-1)*P^k/(I2-f2)*( derivativeA1*I2/(I2-f1)*A2+A1*I2/(I2-f1)*derivativef1/(I2-f1)*A2+A1*I2/(I2-f1)*derivativeA2+derivativeA3+derivativeA4)...
                   + z^k*P^k/(I2-f2)*derivativef2/(I2-f2)*( derivativeA1*I2/(I2-f1)*A2+A1*I2/(I2-f1)*derivativef1/(I2-f1)*A2+A1*I2/(I2-f1)*derivativeA2+derivativeA3+derivativeA4)...
                   + z^k*P^k/(I2-f2)*( secondderivativeA1*I2/(I2-f1)*A2+derivativeA1*I2/(I2-f1)*derivativef1/(I2-f1)*A2+derivativeA1*I2/(I2-f1)*derivativeA2+derivativeA1*I2/(I2-f1)*derivativef1/(I2-f1)*A2+A1*I2/(I2-f1)*derivativef1/(I2-f1)*derivativef1/(I2-f1)*A2+A1*I2/(I2-f1)*secondderivativef1/(I2-f1)*A2+A1*I2/(I2-f1)*derivativef1/(I2-f1)*A2+A1*I2/(I2-f1)*derivativef1/(I2-f1)*derivativef1/(I2-f1)*A2+A1*I2/(I2-f1)*derivativef1/(I2-f1)*derivativeA2 + derivativeA1*I2/(I2-f1)*derivativeA2 + A1*I2/(I2-f1)*derivativef1/(I2-f1)*derivativeA2 + A1*I2/(I2-f1)*secondderivativeA2 +secondderivativeA3+secondderivativeA4);

    
    secondmomentDelay = (pi_I_kron * secondderivativePhiD * onevec)/(pi_I_kron * onevec)+meanDelay;
    
    varDelay = secondmomentDelay-meanDelay^2;
end
