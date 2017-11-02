function [eta] = eigeval(k,l)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
a = 0.05;
b = 0.1;

eta = sqrt((4*(a^l*b+a*b^l)*(2*k+1)*(l+2*k+1)*(l+4*k))...
    /(a*b*(a^(l+1)+b^(l+1))*(l+4*k+2)));

end

