function [thickness] = rho1(a,d,gam,k,l)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
thickness = (8*(d+1)*(2*k+1)*(l+2*k+1)*(l+4*k))./(gam*a*(l+4*k+2))-a;

end

