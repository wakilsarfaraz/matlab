%% This script shows the positivity of the determinant for certain eigenvalues 
% of the RDS in polar coordinates
M = 1000;
rho = 1;
n = 0.3;

alpha = 0.01;
beta = 0.2;
gamma = .1;
d = 10;

k = 0:1:M;



determinant = (gamma.*(beta-alpha)./(beta+alpha)-evaluate(n,k,rho)).*...
    (-gamma.*(beta+alpha).^2-(d+1).*evaluate(n,k,rho))+2*gamma.^2.*beta.*(beta+alpha); 

plot(k,determinant)
xlabel('Eigenvalues')
ylabel('Determinant')
