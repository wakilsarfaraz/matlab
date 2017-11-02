function [thickness] = compthick(d,gam,k,l)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
thickness = (4*(d+1)*(2*k+1)*(l+2*k+1)*(l+4*k)-gam*.5^2*(l+4*k+2))./(gam*.5*(l+4*k+2));

end

