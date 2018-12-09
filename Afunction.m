function [An,derivativeAn,Bn,derivativeBn] = Afunction(z,k,K,T,Anold,derivativeAnold,Bnold,derivativeBnold,P00,P01,P10,P11,Px0,Px1,ncount)

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

    P1x_D = P10.*z.^(k-1)+P11.*z.^(T-1); %The probability of error-free or erroneous NACK   
    
    Px0_C = zeromat;   
    derivativePx0_C = zeromat;
    sum1 = zeromat;
    derivativesum1 = zeromat;
    sum2 = zeromat;
    sum3 = zeromat;
    for j = 0:1000%d-1
       sum1 = sum1 + (z*P1x_D)^j; 
       if j == d-1
           sum2 = sum1;
       end
       if j < T
           sum3 = sum3 + Px1^j;
       end
       derivativesum1 = derivativesum1 + j*z^(j-1)*P1x_D^j;
       Px0_C = Px0_C + z*P01*(z*Px1)^j*Px0;
       derivativePx0_C = derivativePx0_C + P01*(z*Px1)^j*Px0+z*P01*j*z^(j-1)*Px1^j*Px0;
    end  
    
    P01_C = z*P01/(I-z*Px1);
    derivativeP01_C = P01/(I-z*Px1) + z*P01/(I-z*Px1)*Px1/(I-z*Px1);

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    An = z*P00*Anold...
       + z*(z*P01_C*sum1)^ncount * Px0_C;
   
    derivativeAn = P00*Anold + P00*derivativeAnold...
       +( (z*P01_C*sum1)^ncount + ncount*(P01_C*sum1+derivativeP01_C*sum1+P01_C*derivativesum1)*(z*P01_C*sum1)^(ncount-1) ) * Px0_C...
       +z*(z*P01_C*sum1)^ncount * derivativePx0_C;
   
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Bn = P00*Bnold...
        +(P01_C*sum2+P01_C*P1x_D^d/(I-z*Px1^T)*z*sum3)^ncount * Px0_C;
    
    derivativeBn = P00*derivativeBnold...
        +ncount*(P01_C*P1x_D^d/(I-z*Px1^T)*Px1^T/(I-z*Px1^T)*z*sum3+P01_C*P1x_D^d/(I-z*Px1^T)*sum3)*(P01_C*sum2+P01_C*P1x_D^d/(I-z*Px1^T)*z*sum3)^(ncount-1) * Px0_C;

end