function [P00,P01,P10,P11,P0x,Px0,P1x,Px1,P_kron,pi_kron,pi_I_kron] = TransitionProbabilityMatrices(eps_Gf,eps_Bf,epsf,eps_Gr,eps_Br,epsr,rf,rr,NACK,K)
%Computation of the transition probability matrices                                                 
    if K == 2      
        if rf ~= 0 
            e_vecf = [eps_Gf, eps_Bf];         %error vector for forward channel
            e_vecr = [eps_Gr, eps_Br];         %error vector for reverse channel

            qf = rf*(eps_Gf-epsf)/(epsf-eps_Bf);  %transition probability from G to B  
            qr = rr*(eps_Gr-epsr)/(epsr-eps_Br);

            %P = [GG, GB; BG, BB]: r is the transition probability from B to G
            Pf = [1-qf, qf; rf, 1-rf];           %state transition matrix for forward channel
            Pr = [1-qr, qr; rr, 1-rr];           %state transition matrix for reverse channel
            %1/r represents the average error burst

            pivecf = [rf/(rf+qf), 1-rf/(rf+qf)];   %stationary distribution of the MC
            pivecr = [rr/(rr+qr), 1-rr/(rr+qr)];

        else %rf == 0   %MARKOVian errors, always in the bad state
            e_vecf = epsf;
            e_vecr = epsr;

            qf = 0;
            qr = 0;

            Pf = 1;
            Pr = 1;

            pivecf = 1;
            pivecr = 1;

        end
        
    else
        e_vecf = epsf;
        e_vecr = epsr;

        qf = 0;
        qr = 0;

        Pf = 1;
        Pr = 1;

        pivecf = 1;
        pivecr = 1;
    end

    %eps= pivec*e_vec';        

    %P = P0 + P1
    P1f = Pf * diag(e_vecf);       %error probability matrix of HMM
    P0f = Pf * diag(1-e_vecf);     %success probability matrix of HMM
    
    P1r = Pr * diag(e_vecr);
    P0r = Pr * diag(1-e_vecr);

    %eps = pivec*P1*ones(K,1);
    
    %P1f = P1;       
    %P0f = P0;
    %P1r = P1f;
    %P0r = P0f;    
        
    %Kronecker product of two matrices: Pfr = kron(Pf,Pr)
    
    if NACK == 1   %with NACK (regular scheme)
        P00 = kron(P0f,P0r);    %error-free ACK
        P01 = kron(P0f,P1r);    %erroneous ACK
        P10 = kron(P1f,P0r);    %error-free NACK
        P11 = kron(P1f,P1r);    %erroneous NACK
    else           %without NACK (simplified scheme)
        P00 = kron(P,P0r);      %error-free ACK
        P01 = kron(P,P1r);      %erroneous ACK
        P10 = zeros(2*K,2*K);   %error-free NACK
        P11 = zeros(2*K,2*K);   %erroneous NACK
    end
            
    P0x = P00+P01;    %probability matrix of success in the forward channel
    Px0 = P00+P10;    %probability matrix of success in the reverse channel
    P1x = P10+P11;    %probability matrix of error in the forward channel
    Px1 = P01+P11;    %probability matrix of error in the reverse channel
    
    %Kronecker product of forward and reverse channels "Pf \otimes Pr"
    
    P_kron = P00 + P01 + P10 + P11;
    
    pi_kron = kron(pivecf,pivecr);

    pi_If = pivecf*P0f; 
    pi_Ir = pivecr*P0r; 
    
    pi_I_kron = kron(pi_If,pi_Ir);  
    
    
end
