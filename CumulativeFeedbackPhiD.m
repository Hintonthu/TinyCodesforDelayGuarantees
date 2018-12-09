function [PhiD, phiD, meanDelay, Throughput] = CumulativeFeedbackPhiD(z,k,K,T,P_kron,pi_I_kron,P00,P01,P10,P11,Px0,Px1,M,N,symbolic)
%Matrix-generating function of delay, coding

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
    
   
    
    RTT = k+M-1;   %k = RTT-M+1; %time needed to get the first feedback after the last packet is transmitted
    d = T-RTT;     %time to timeout    

   
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
    I2 = I; %[I, zeromat;zeromat, I];

    
    A = I-z^k*P10*P_kron^(k-1)-z^T*P11*P_kron^(T-1);
    B = I-z*P10*P_kron^(k-1)-z*P11*P_kron^(t-1);
    C = I-Px1;
    
    %b*inv(A) = b/A         %inv(A)*b = A\b
    
    P = P_kron;
    
    % ANALYTICAL EXPRESSION FOR DELAY   
    P10C = P10*P10+P10*P01+P01*P10;
    P11C = P11*P11+P11*P10+P10*P11;
    P00C = P00*P00;
    P01C = P01*P01+P01*P00+P00*P01;
    P0 = P00*P10+P10*P00...
       + P01*P00+P00*P01...+ P11*P01+P01*P11... %
       + P00*P11+P11*P00;
      
    Px0C = P00C+P10C; %Px0*Px0;
    Px1C = P01C+P11C; %Px1*Px1;
                          
    
    P1x_D_1 = z^(T-d)*P^(T-d)*(z*P10C+z*P11C*z^d*P^d);
    P1x_D_2 = z^(T-d-1)*P^(T-d-1)*(z*P10+z*P11*z^(d+1)*P^(d+1));
    PhiD = z^k*P^k/(I2-P1x_D_1)*(z*P00C+z*P01C/(I2-z*Px1C)*z*Px0C+z*P0/(I2-P1x_D_2)*(z*P00+z*P01/(I2-z*Px1)*z*Px0));
      
    phiD = (pi_I_kron * PhiD * onevec)/(pi_I_kron * onevec);
    
    %Derivative of phiD: evaluate at z=1 to obtain mean delay
    
    derivativeP1x_D_1 = (T-d)*z^(T-d-1)*P^(T-d)*(z*P10C+z*P11C*z^d*P^d)...
                      + z^(T-d)*P^(T-d)*(P10C+(d+1)*P11C*z^d*P^d);
    derivativeP1x_D_2 = (T-d-1)*z^(T-d-2)*P^(T-d-1)*(z*P10+z*P11*z^(d+1)*P^(d+1))...
                      + z^(T-d-1)*P^(T-d-1)*(P10+(d+1)*P11*z^d*P^(d+1));
                    
    derivativePhiD = k*z^(k-1)*P^k/(I2-P1x_D_1)*(z*P00C+z*P01C/(I2-z*Px1C)*z*Px0C+z*P0/(I2-P1x_D_2)*(z*P00+z*P01/(I2-z*Px1)*z*Px0))...
                   + z^k*P^k/(I2-P1x_D_1)*derivativeP1x_D_1/(I2-P1x_D_1)*(z*P00C+z*P01C/(I2-z*Px1C)*z*Px0C+z*P0/(I2-P1x_D_2)*(z*P00+z*P01/(I2-z*Px1)*z*Px0))...
                   + z^k*P^k/(I2-P1x_D_1)...
                   *(P00C+2*z*P01C/(I2-z*Px1C)*Px0C+z*P01C/(I2-z*Px1C)*Px1C/(I2-z*Px1C)*z*Px0C+P0/(I2-P1x_D_2)*(z*P00+z*P01/(I2-z*Px1)*z*Px0)...
                   +z*P0/(I2-P1x_D_2)*derivativeP1x_D_2/(I2-P1x_D_2)*(z*P00+z*P01/(I2-z*Px1)*z*Px0)...
                   +z*P0/(I2-P1x_D_2)*(P00+2*z*P01/(I2-z*Px1)*Px0+z*P01/(I2-z*Px1)*Px1/(I2-z*Px1)*z*Px0));

    %Mean delay per M packets (bucket)           
    meanDelay = (pi_I_kron * derivativePhiD * onevec)/(pi_I_kron * onevec); 
    
    
    % ANALYTICAL EXPRESSION FOR THROUGHPUT 
    P1x_T_1 = z*(P10C+P11C*P^d*P^d)*P^(T-d);
    P1x_T_2 = z*(P10+P11*P^(d-1))*P^(T-d-1);
    if symbolic
        sum1 = (1-Px1C^T)/(1-Px1C);
        sum2 = P01C*(1-Px1C^d)/(1-Px1C)*Px0C;
        sum3 = P01*(1-Px1^(d-1))/(1-Px1)*(Px0);
        sum4 = (1-Px1^T)/(1-Px1);
    else
        sum1 = I2;
        sum2 = zeromat;
        sum3 = zeromat;
        sum4 = I2; %this was zeromat before
        for i = 1:T-1
            sum1 = sum1 + Px1C^i; %1st layer where 2DoFs required      
            if i<=d
                sum2 = sum2 + P01C*Px1C^(i-1)*Px0C;     %1st layer where 2DoFs required      
            end
            if i<=d-1
                sum3 = sum3 + (P01*Px1^(i-1)*Px0);      %2nd layer where 1DoF required  
            end
            sum4 = sum4 + Px1^i;  %2nd layer where 1DoF required    
        end
    end
    Phitau = z*P^(T-d-1)/(I2-P1x_T_1)*(P01C*Px1C^d/(I2-z*Px1C^T)*z*sum1*Px0C...
                                     +P00C+sum2...
                                     +P0/(I2-P1x_T_2)*(P00+sum3+P01*Px1^(d-1)/(I2-z*Px1^T)*z*sum4*Px0));    
       
                                 
    derivativeP1x_T_1 = (P10C+P11C*P^d*P^d)*P^(T-d);
    derivativeP1x_T_2 = (P10+P11*P^(d-1))*P^(T-d-1);
                                 
    derivativePhitau = P^(T-d-1)/(I2-P1x_T_1)*(P01C*Px1C^d/(I2-z*Px1C^T)*z*sum1*Px0C...
                                     +P00C+sum2...
                                     +P0/(I2-P1x_T_2)*(P00+sum3+P01*Px1^(d-1)/(I2-z*Px1^T)*z*sum4*Px0))...
                     + z*P^(T-d-1)/(I2-P1x_T_1)*derivativeP1x_T_1/(I2-P1x_T_1)*(P01C*Px1C^d/(I2-z*Px1C^T)*z*sum1*Px0C...
                                     +P00C+sum2...
                                     +P0/(I2-P1x_T_2)*(P00+sum3+P01*Px1^(d-1)/(I2-z*Px1^T)*z*sum4*Px0))...
                     + z*P^(T-d-1)/(I2-P1x_T_1)...
                     *(P01C*Px1C^d/(I2-z*Px1C^T)*(Px1C^T)/(I2-z*Px1C^T)*z*sum1*Px0C+P01C*Px1C^d/(I2-z*Px1C^T)*sum1*Px0C...                                  
                                     +P0/(I2-P1x_T_2)*derivativeP1x_T_2/(I2-P1x_T_2)*(P00+sum3+P01*Px1^(d-1)/(I2-z*Px1^T)*z*sum4*Px0)...
                                     +P0/(I2-P1x_T_2)*(P01*Px1^(d-1)/(I2-z*Px1^T)*Px1^T/(I2-z*Px1^T)*z*sum4*Px0+P01*Px1^(d-1)/(I2-z*Px1^T)*sum4*Px0));                 
                                 
     
    %Total transmission time for M packets                             
    transmissionTime = (pi_I_kron * derivativePhitau * onevec)/(pi_I_kron * onevec); 
    
    Throughput = 1/transmissionTime;    

end
%taylor(meanDelay,e) %5th order Taylor approximation
%taylor(Throughput,e)
%polyfit(e,Throughput,3)
