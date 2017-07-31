%% Coefficient of the eigenvalues
lmax = 10;
N = 100;
a = 1.5;
b = 3;

l = linspace(-lmax,lmax,N);
f = zeros(1,N);
for k = 1 : N
f(k)=(a.^l(k)*b+a*b.^l(k))/(a*b*(a.^(l(k)+1)+b.^(l(k)+1)));
end

% plot(l,f,'-','linewidth',2)
% legend('a=0.5, b=1.0','a=0.6,b=1.2','a=0.7,b=1.4','a=0.8,b=1.6','a=0.9,b=1.8','a=1.0,b=2.0','a=1.5,b=3.0')
% title('The variation of the domain depedent factor of \eta^2 for different \rho','fontsize',10)
% xlabel('l','fontsize',18)
% ylabel('\eta^2','fontsize',18)
% hold on


%    set(findobj('type','legend'),'fontsize',18)
%    set(findobj('type','axes'),'fontsize',16)