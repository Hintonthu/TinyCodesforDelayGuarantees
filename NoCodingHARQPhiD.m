function [PhiD, phiD, meanDelay, Throughput] = NoCodingHARQPhiD(z,k,K,T,eps,r,NACK,symbolic)                                                    
%Matrix-generating function of delay, no coding HARQ

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
    
    A = I-z^k*P10*P_kron^(k-1)-z^T*P11*P_kron^(T-1);
    B = I-P10*P_kron^(k-1)-P11*P_kron^(t-1);
    %C = I-Px1;
    
    
    %b*inv(A) = b/A         %inv(A)*b = A\b
    
    
    % ANALYTICAL EXPRESSION FOR DELAY

    PhiD = z^(k-1)*P_kron^(k-1)/A.*(z*P00+z^2*(P01/(I-z*Px1)) *Px0);
    %PhiD = z^(k-1)*P_kron^(k-1)*C1.*(z*P00+z^2*P01*C2); %updated
    
             
    phiD = pi_I_kron * PhiD * ones(2*K,1)/(pi_I_kron * onevec);
    
    %Derivative of phiD: evaluate at z=1 to obtain mean delay
                       
    derivativePhiD = (P_kron^(k-1)/B)...
                     *(k*P10*P_kron^(k-1)+t*P11*P_kron^(t-1))...
                     /B...                          %derivative of 1st term
                     *(P00+P01*D1)...               %2nd term
                     +(P_kron^(k-1)/B)*...          %1st term
                     (k*P00+(k+1)*P01*D1...         %derivative of 2nd term
                     +P01*D2);         
    %OR THE LAST TERM OF THE LAST LINE: P01*D1 * D3^(1/cmax) * D1   
    %derivativePhiD = P_kron^(k-1)*C1derivative.*(z^k*P00+z^(k+1)*P01*C2)...  
    %               + P_kron^(k-1)*C1.*(k*z^(k-1)*P00+(k+1)*z^k*P01*C2+z^(k+1)*P01*C2derivative); %updated
                    
    meanDelay = (pi_I_kron * derivativePhiD * onevec)/(pi_I_kron * onevec);
    
    
    % ANALYTICAL EXPRESSION FOR THROUGHPUT

    Phitau = z*P_kron^(k-1)/(I-z*P10*P_kron^(k-1)-z*P11*P_kron^(T-1))...
           *(P00+P01*sum1 ...
           +P01*sum2*z*sum3);

    derivativePhitau = (P_kron^(k-1)/(I-P10*P_kron^(k-1)-P11*P_kron^(T-1))...
        +P_kron^(k-1)/(I-P10*P_kron^(k-1)-P11*P_kron^(T-1))*(P10*P_kron^(k-1)+P11*P_kron^(T-1))/(I-P10*P_kron^(k-1)-P11*P_kron^(T-1)) )...
       *(P00+P01*sum1...
       +P01*sum2*z*sum3)...
       +P_kron^(k-1)/(I-P10*P_kron^(k-1)-P11*P_kron^(T-1))...
       *(P01*sum2derivative*sum3+P01*sum2*sum3);    

    

   
    transmissionTime = (pi_I_kron * derivativePhitau * onevec)/(pi_I_kron * onevec);

    
    Throughput = 1/transmissionTime;   
              
end