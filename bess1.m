%% Coefficient of the eigenvalues with parameter introduced
rhomax = 10;
N = 100;
a = 0.5;
l = 0.3;

rho = linspace(-rhomax,rhomax,N);
f = zeros(1,N);
for k = 1 : N
f(k)=(a^(l-1)*(1+a*(rho(k)+a)^(l-1)))/(a^(l+1)+(rho(k)+a)^(l+1));
end

 plot(l,f,'-','linewidth',2)
% legend('a=0.5, b=1.0','a=0.6,b=1.2','a=0.7,b=1.4','a=0.8,b=1.6','a=0.9,b=1.8','a=1.0,b=2.0','a=1.5,b=3.0')
% title('The variation of the domain depedent factor of \eta^2 for different \rho','fontsize',10)
% xlabel('l','fontsize',18)
% ylabel('\eta^2','fontsize',18)
% hold on


%    set(findobj('type','legend'),'fontsize',18)
%    set(findobj('type','axes'),'fontsize',16)