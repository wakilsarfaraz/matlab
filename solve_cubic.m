function [] = solve_cubic()
w=10:2:30;
a=4*((9*w)-7);
b=(-2)*((18*w)-16).*(w+1);
c=3*((w+1).^2).*((3*w)-1);
d=(-3)*((1+w).^3);
f = zeros(1,length(w));
for ii = 1:numel(w)
    g = roots([a(ii),b(ii),c(ii),d(ii)]);
    f(ii) = g(g<1);
end
ptr=2./(1+w);
u=(1+w)./2;
l=(u.*(1-ptr.*f).^3)./3;
plot(w,l);