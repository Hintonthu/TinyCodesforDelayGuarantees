%Gaussian Tail
clear all; clc; %close all;

T = 15;

C = 1;%constant
v = 3*10^-4; %constant

%P(|X|>d) <= C e^{-vd^2}
CCDF_Gaussian = @(d) C*exp(-v*d.^2);

n_min = 1;
n_max = 5*T;
n_set = linspace(n_min,n_max,20);
figure    
plot(n_set,CCDF_Gaussian(n_set),'k','linewidth',2) 
xlab = 'Value of delay, d'; 
ylab = 'CCDF, $P(D_{\sf CF-ARQ}>d$)';
box on;     set(gca,'FontSize',20) 
xlhand = get(gca,'xlabel'); xlabel(xlab,'Interpreter','latex'); set(xlhand,'fontsize',20) 
ylhand = get(gca,'ylabel'); ylabel(ylab,'Interpreter','latex'); set(ylhand,'fontsize',20)
xaxis = T;    yaxis = 0.1;
text(xaxis,yaxis,['$C$=' num2str(C) ', $v$=' num2str(v)],'fontsize',20,'Interpreter','latex');
