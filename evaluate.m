function [eta_squared] = evaluate(n,k,rho)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
eta_squared = - (4.*(2.*k+1).*(n+2.*k+1).*(n+k))./(rho.^2.*(n+4.*k+2));

end

