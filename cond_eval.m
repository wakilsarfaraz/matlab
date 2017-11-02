function [rho_cond] = cond_eval(n,k,d,gam)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
rho_cond =2* sqrt(((d+1)*(2*k+1)*(n+2*k+1)*(n+4*k))/(gam*(n+4*k+2)));

end

