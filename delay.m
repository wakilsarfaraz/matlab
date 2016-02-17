%%This script solve a system of coupled reaction diffusion equations with
%%time delay term, and zero Neumann B.C.
% Wakil Sarfaraz 17/02/2016


%clear all; 
close all; 
clc;

xmax = 1;
tm = 1;
dt = 0.1;
M = tm/dt;

Dh = 0.01;
Dp = 0.025;

a1 = 1.2;
a2 = 0.3;
b = 0.2;
c1 = 0.9;
c2 = 0.42;
k1 = 1;
k2 = 2.1;

tau = 0.5; % Delay Parameter

N = 10;
X = linspace(0,xmax,N+1);
T = linspace(0,tm,M+1);

[x,y] = meshgrid(X,X);

x = x(:);
y = y(:);

NNODES = (N+1)^2;
H = zeros(NNODES,1);
P = zeros(NNODES,1);

NTRI = 2*N^2;  
LNODES = zeros(NTRI,3); 
for i = 1:N
    for j = 1: N
        LNODES(i+2*(j-1)*N,1) = i+(j-1)*(N+1);
        LNODES(i+2*(j-1)*N,2) = i+j*(N+1); 
        LNODES(i+2*(j-1)*N,3) = (i+1)+(j-1)*(N+1); 
        
        LNODES(i+N+2*(j-1)*N,1) = i+1+j*(N+1); 
        LNODES(i+N+2*(j-1)*N,2) = (i+1)+(j-1)*(N+1);
        LNODES(i+N+2*(j-1)*N,3) = i+j*(N+1);
                                   
    end
end


GH = zeros(NNODES,2);
GP = zeros(NNODES,2);

for i = 1 : NNODES
    H(i) = H(i)+0.0001*pi^2*sin(10*xmax*pi*sin(5*pi*x(i))+10*xmax*pi*sin(5*pi*y(i)));%*sin(0.5*pi*x(i));
    P(i) = P(i)+0.0001*pi^2*sin(15*pi*x(i))*cos(15*pi*y(i));
end

for i = 1 : NNODES
    if (y(i)==0 && x(i)>=0 && x(i) <= xmax && i<=N)
        GH(i,2)= (H(LNODES(i,2))-H(LNODES(i,1)))/(y(LNODES(i,2))-y(LNODES(i,1)));
    else

    end
end
trisurf(LNODES,x,y,P(:,:))

GH











































































