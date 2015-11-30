%%System of nonlinear odes
%Author Wakil Sarfaraz 17/11/2015

a1=1.2; a2=0.3; b=0.2; c1=0.9; c2=0.42; k1=1;


f = @(t,x) [x(1)*a1-b*(x(1))^2-(c1*x(1)*x(2)); a2*x(2)-(c2*(x(2))^2)/(k1+x(1))];

[t,xa] = ode45(f,[0 30],[0.5 0.5]);


plot(t,xa(:,1),'-r',t,xa(:,2),'-.b')
grid on
title('Absence of diffusion with parameters of Sec3')
xlabel('t'), ylabel('H(t) & P(t)')
legend('Prey denoted by H(t)','Preditor denoted by P(t)')


