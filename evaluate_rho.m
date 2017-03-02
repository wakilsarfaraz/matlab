function [rho] = evaluate_rho(d,g,n,k)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
rho= 2*sqrt(((d+1)*(2*k+1)*(n+2*k+1)*(n+4*k))/(g*(n+4*k+2)));

end