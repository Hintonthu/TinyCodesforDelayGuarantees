function [inverse1,inverse2] = MatrixInverse(A,z)

    s = size(A,1);
    inverse1 = eye(s,s)/(eye(s,s)-z*A)*A/(eye(s,s)-z*A);
    
    
    inverse2 = zeros(s,s); 
    for j = 1:100
       inverse2 = inverse2 + j*z^(j-1)*A^j; 
    end

end