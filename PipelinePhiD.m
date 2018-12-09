function [PhiD, phiD, meanDelay, Throughput] = PipelinePhiD(z,k,K,T,P_kron,pi_I_kron,P00,P01,P10,P11,Px0,Px1,M,N)
%Matrix-generating function of delay, coding

    %M packets transmitted
    %n DoFs required

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
    C = I-Px1;
    
    %b*inv(A) = b/A         %inv(A)*b = A\b
    
    
    % ANALYTICAL EXPRESSION FOR DELAY   
    [An,derivativeAn,Bn,derivativeBn] = Afunction(z,k,K,T,1,0,1,0,P00,P01,P10,P11,Px0,Px1,1);  %A0=1, %B0=1
    for ncount = 2:N
        [An,derivativeAn,Bn,derivativeBn] = Afunction(z,k,K,T,An,derivativeAn,Bn,derivativeBn,P00,P01,P10,P11,Px0,Px1,ncount);
    end
    
    PhiD = z^(k-1)*P_kron^(k-1)/A * An;
      
    phiD = (pi_I_kron * PhiD * onevec)/(pi_I_kron * onevec);
    
    %Derivative of phiD: evaluate at z=1 to obtain mean delay
                    
    derivativePhiD = (k-1)*P_kron^(k-1)/A * An...
                   + P_kron^(k-1)/A*(k*P10*P_kron^(k-1)+T*P11*P_kron^(T-1))/A* An...
                   + P_kron^(k-1)/A* derivativeAn;

    meanDelay = (pi_I_kron * derivativePhiD * onevec)/(pi_I_kron * onevec);
    meanDelay = meanDelay/M; %Mean delay per packet
    
    
    % ANALYTICAL EXPRESSION FOR THROUGHPUT                           
    Phitau = z*P_kron^(k-1)/B * Bn;    
       
    derivativePhitau = P_kron^(k-1)/B * Bn...
                     + P_kron^(k-1)/B * (P10*P_kron^(k-1)+P11*P_kron^(t-1))/B * Bn...
                     + P_kron^(k-1)/B * derivativeBn;
        
    transmissionTime = (pi_I_kron * derivativePhitau * onevec)/(pi_I_kron * onevec); %Transmission time for M packets
    
    Throughput = 1/(transmissionTime/N);   %/M

end
