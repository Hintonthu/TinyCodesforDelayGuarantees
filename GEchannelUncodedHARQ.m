
%Uncoded HARQ, evolution of GE channel (single realization)

function [ChannelStatef,ChannelStater,Pkron] = GEchannelUncodedHARQ(q,q_HARQ,r,total_no)

    %assume m (number of attempts) ranges from 1:total_no

    m_set= 1:total_no; 
    q_sequence = q_HARQ(m_set); 

    for m = 1:m_set
        q_HARQ = q_sequence(m); %instead of having fixed eps for all channel uses, we have decreasing eps
        
        P = [1-q, q; r, 1-r];
        Pkron = kron(P,P);
        pi_G_f = r/(r+q); %probability of being in State G
        pi_B_f = q/(r+q); %probability of being in State B
        
        pi_G_r = r/(r+q_HARQ); %probability of being in State G
        pi_B_r = q_HARQ/(r+q_HARQ); %probability of being in State B
        
        %eps = pi_B*eps_B;
        %eps_B = 1/(1+r/q)/pi_B;
       
        %pi_I_G = pi_G*(1-q)+pi_B*r;
        %pi_I_B = (pi_G*q+pi_B*(1-r))*(1-epsB);
   
        ChannelStatef = zeros(1,total_no);
        ChannelStater = zeros(1,total_no);
    
        goodf = rand(1) > pi_B_f;
        goodr = rand(1) > pi_B_r;          

        for size = 1:total_no

            if goodf == 1 && goodr == 1        
                goodf = rand(1) > q ; 
                goodr = rand(1) > q_HARQ ; 

            elseif goodf == 1 && goodr == 0
                goodf = rand(1) > q ; 
                goodr = rand(1) > 1-r; 

            elseif goodf == 0 && goodr == 1
                goodf = rand(1) > 1-r; 
                goodr = rand(1) > q_HARQ; 

            elseif goodf == 0 && goodr == 0
                goodf = rand(1) > 1-r; 
                goodr = rand(1) > 1-r; 

            end
            ChannelStatef(size) = goodf;
            ChannelStater(size) = goodr;
        end
        
    end
    
                
end
        
