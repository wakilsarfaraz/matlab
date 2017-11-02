%%Eigenvalues Matrix
k = 12;
nl = 30;
ldecim = -11.3;
val = ones(1,nl);

for i = 1 : nl
    val(i) = (ldecim+i-1)*val(i);
end
M = zeros(k,nl);
for i = 1 : k
    for j = 1 : nl
        M(i,j) = M(i,j)+eigeval(i,val(j));
    end
end
for i = 1 : k
    plot(val,M(i,:),'*--','linewidth',2)
    axis tight
    hold on
    pause(1e-10)
end
xlabel('Order of Bessel`s equation l','fontsize',18)
ylabel('Eigenvalue \eta_{k,l}','fontsize',18)
% xlim([-11.3 11.3])
% ylim([-50 90])
legend('k = 1','k = 2','k = 3','k = 4','k = 5','k = 6',...
    'k = 7','k = 8','k = 9','k = 10','k = 11','k = 12','location','NorthEast')
title('Variation of \eta_{k,l} with a = 0.9 and b = 1','fontsize',20)
   set(findobj('type','legend'),'fontsize',16)
 set(findobj('type','axes'),'fontsize',20)
