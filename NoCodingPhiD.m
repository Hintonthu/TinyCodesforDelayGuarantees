function [PhiD, phiD, meanDelay,varDelay, Throughput] = NoCodingPhiD(z,k,K,T,P_kron,pi_I_kron,P00,P01,P10,P11,Px0,Px1,symbolic)
%Matrix-generating function of delay, no coding

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
        
    A = I-z^k*P10*P_kron^(k-1)-z^T*P11*P_kron^(T-1);
    B = I-z*P10*P_kron^(k-1)-z*P11*P_kron^(t-1);
    Btau = I-z*P10*P_kron^(k-1)-z*P11*P_kron^(T-1);
    C = I-Px1;
    
    %b*inv(A) = b/A         %inv(A)*b = A\b

    
    % ANALYTICAL EXPRESSION FOR DELAY
    %P1x_D = @(z) P10.*z.^(k-1)+P11.*z.^(T-1); %The probability of error-free or erroneous NACK   

   
    PhiD = z^(k-1)*P_kron^(k-1)/A.*(z*P00+z^2*(P01/(I-z*Px1))*Px0);    
        
    %Derivative of phiD: evaluate at z=1 to obtain mean delay                       
    derivativePhiD = (P_kron^(k-1)/B)...
                 *(k*P10*P_kron^(k-1)+t*P11*P_kron^(t-1))...
                 /B...                          %derivative of 1st term
                 *(P00+P01/C*Px0)...            %2nd term
                 +(P_kron^(k-1)/B)*...          %1st term
                 (k*P00+(k+1)*P01/C*Px0...      %derivative of 2nd term
                  +P01/C*Px1/C*Px0);
 
    
    phiD = pi_I_kron * PhiD * onevec/(pi_I_kron * onevec);
    

    meanDelay = (pi_I_kron * derivativePhiD * onevec)/(pi_I_kron * onevec);
    
    
    % ANALYTICAL EXPRESSION FOR THROUGHPUT
    if symbolic
        sum1 = (1-Px1^d)/(1-Px1);
        sum2 = (1-Px1^T)/(1-Px1);       
     else
        sum1 = zeromat;
        sum2 = zeromat;
        for i = 0:T-1
            if i <= d-1
                sum1 = sum1 + Px1^i;
            end
            sum2 = sum2 + Px1^i;
        end   
    end
                            
    Phitau = z*P_kron^(k-1)/Btau...
           *(P00 + P01*sum1*Px0 + P01*Px1^d/(I-z*Px1^T) * z * sum2 * Px0);
            
    derivativePhitau = (P_kron^(k-1)/Btau+P_kron^(k-1)/Btau*(P10*P_kron^(k-1)+P11*P_kron^(T-1))/Btau)...
           *(P00 + P01*sum1*Px0 + P01*Px1^d/(I-Px1^T) * sum2 * Px0)...
           +P_kron^(k-1)/Btau...
           *(P01*Px1^d/(I-Px1^T)*Px1^T/(I-Px1^T) * sum2 * Px0 + P01*Px1^d/(I-Px1^T) * sum2 * Px0);
            
    %OLD THROUGHPUT DERIVATIVE FUNCTION                                 
%     derivativePhitau = (P_kron^(k-1)/Btau...
%            +P_kron^(k-1)/Btau*(P10*P_kron^(k-1)+P11*P_kron^(T-1))/Btau)...
%            *(P00+P01*sum1*Px0+P01*Px1^d/(I-Px1^T) * sum2 * Px0)...
%            +P_kron^(k-1)/Btau...
%            *(P01*Px1^d/(I-Px1^T)*Px1^T/(I-Px1^T) * sum2 * Px0 + P01*Px1^d/(I-Px1^T) * sum2 * Px0);

       
    transmissionTime = (pi_I_kron * derivativePhitau * onevec)/(pi_I_kron * onevec);
    
    Throughput = 1/transmissionTime;   
    
    
    %VARIABILITY OF DELAY: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    A = P10*P_kron^(k-1)+P11*P_kron^(T-1);
    Afirstder = k*P10*P_kron^(k-1)+T*P11*P_kron^(T-1);
    Asecondder = k*(k-1)*P10*P_kron^(k-1)+T*(T-1)*P11*P_kron^(T-1);
    
    B = ((I-A)\Afirstder)/(I-A);
    
    Bfirstder = 2*((I-A)\Afirstder)*((I-A)\Afirstder)/(I-A) + ((I-A)\Asecondder)/(I-A);
    
    C = ((I-Px1)\Px1)/(I-Px1);
    
    Cfirstder = 2*((I-Px1)\Px1)*((I-Px1)\Px1)/(I-Px1);
    
    %b*inv(A) = b/A         %inv(A)*b = A\b
    
    %Second derivative of phiD: evaluate at z=1 to obtain the second moment of delay
    secondderivativePhiD = ((k-1)*(k-2)*P_kron^(k-1)/(I-A)+(k-1)*P_kron^(k-1)*B+(k-1)*P_kron^(k-1)*B+P_kron^(k-1)*Bfirstder)...
                         * (P00+P01*((I-Px1)\Px0))...
                         + 2*((k-1)*P_kron^(k-1)/(I-A)+P_kron^(k-1)*B)...
                         * (P00+2*(P01/(I-Px1))*Px0+P01*C*Px0)...                      
                         + P_kron^(k-1)/(I-A)...
                         *(2*(P01/(I-Px1))*Px0+4*P01*C*Px0+P01*Cfirstder*Px0);
    
    %Second Moment of delay
    secondmomentDelay = (pi_I_kron * secondderivativePhiD * onevec)/(pi_I_kron * onevec)+meanDelay;
    
    varDelay = secondmomentDelay-meanDelay^2;
end