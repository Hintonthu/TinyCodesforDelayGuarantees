 clear all; clc; close all

 alphaCF = @(x)(1+3*x-2*x^2+20*x^3-18*x^4+28*x^5-60*x^6+72*x^7-40*x^8+8*x^9)/(2*(2-x)*(1-x+4*x^2-2*x^3)*(1+x-2*x^2+2*x^3));

 betaCF = @(x,T,k)(1-x)^2*(1-2*x+4*x^2)/(2*(2-x)*(1-x+4*x^2-2*x^3)*(x^k*(2-x)^k*(1-2*x+2*x^2)^(k+1)...
                         -x^(-(T-k))*(2-x)^(-(T-k))*(1-2*x+2*x^2)^(-(T-k-1))));

 
 xset = linspace(0,0.5,100);
 
 T = 100;
 k = 10;
 
 for i = 1:100
     x= xset(i);
     
    nom(i) = (1+x-2*x^2+2*x^3)/(2-x); 
     
     
    c(i) = alphaCF(x)-betaCF(x,T,k) ;
    
    
    e(i) = x^(T-k)*(1-x)/(1-x^T);
    
    thputgap(i) = (1+x-2*x^2+2*x^3)/(2-x)/(c(i)+x^(T-k)*(1-x)/(1-x^T)) - (1-x)/(1+x^(T-k+1)*(1-x)/(1-x^T));
 end
 figure
 plot(xset,c);
 
% figure
% plot(xset,1./nom)

%figure
%plot(xset,1./(c.*xset))


figure
plot(xset,c./nom.*(1-xset))

% 
% figure
% plot(xset,x.*c)
% 
% figure
% plot(xset,thputgap)