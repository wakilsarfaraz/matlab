%% Coefficient of the eigenvalues with parameter introduced


rhomax = 10;
N = 100;
a = 2;
l = 2.3;

rho = linspace(0,rhomax,N);
f = zeros(1,N);
for k = 1 : N
f(k)=(a^(l-1)*(1+a*(rho(k)+a)^(l-1)))/(a^(l+1)+(rho(k)+a)^(l+1));
end

 plot(rho,f,'-','linewidth',2)
legend('a=0.5','a=1.0','a=1.5','a=2.0')%,'','','')
title('f(\rho) for different values of a')%,'fontsize',14)
xlabel('\rho=b-a','fontsize',18)
ylabel('f(\rho)','fontsize',18)
 hold on


   set(findobj('type','legend'),'fontsize',18)
   set(findobj('type','axes'),'fontsize',18)
%    
