function [shell_eigen ] = eta(a,rho,k,l)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
shell_eigen = sqrt((a.^(l-1)+(rho+a).^(l-1).*(2*k+1)*(l+2*k+1)*(l+4*k))./...
    (a.^(l+1)+(rho+a).^(l+1).*(l+4*k+2)));

end

