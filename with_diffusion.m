%% This script solves a Homogeneous Heat equation in [0 1]X[0 1] with zero Dirichlet B.C.
% Wakil Sarfaraz   30/11/15
clear all;
tic
xmax = 1;
tm = 1;
dt = 0.0001;
M = tm/dt;
N = 50; 
X = linspace(0,xmax,N+1);
T = linspace(0, tm,M+1); 
[x, y] = meshgrid(X,X); 
x =x(:);  
y =y(:); 
NNODES = (N+1)^2;
Q = ones(NNODES,1);
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

SP = sparse(NNODES, NNODES);
SPM = sparse(NNODES, NNODES);
LV = zeros(NNODES,1);         

for n = 1: NTRI
    r1 = [x(LNODES(n,1)) y(LNODES(n,1))];
    r2 = [x(LNODES(n,2)) y(LNODES(n,2))];
    r3 = [x(LNODES(n,3)) y(LNODES(n,3))];
    J = [r2(1)-r1(1) r2(2)-r1(2); r3(1)-r1(1) r3(2)-r1(2)];

Astiff = (1/(2*det(J)))* [(r2-r3)*(r2-r3)' (r2-r3)*(r3-r1)' (r2-r3)*(r1-r2)';... 
           (r2-r3)*(r3-r1)' (r3-r1)*(r3-r1)' (r3-r1)*(r1-r2)';...
           (r2-r3)*(r1-r2)' (r3-r1)*(r1-r2)' (r1-r2)*(r1-r2)'];

   
 Amass = det(J)*[1/12 1/24 1/24; 1/24 1/12 1/24; 1/24 1/24 1/12]; 
 
       for i = 1 : 3 
           for j = 1:3 
               SP(LNODES(n,i),LNODES(n,j)) =SP(LNODES(n,i),LNODES(n,j))+ Astiff(i,j);
               SPM(LNODES(n,i),LNODES(n,j)) =SPM(LNODES(n,i),LNODES(n,j))+ Amass(i,j);
           end
       end      
end
TMatrix =  (dt*SP+SPM);
for i = 1: NNODES
    if (x(i)==xmax && y(i)>=0 && y(i) <= xmax)  
        LV(i) = 0;
        TMatrix(i,:) = 0;
        SPM(i,:) = 0;
        TMatrix(i,i) =1 ;
 
    elseif (x(i)==0 && y(i)>=0 && y(i) <= xmax)  
        LV(i) = 0;
        TMatrix(i,:) = 0;
        SPM(i,:) = 0;
        TMatrix(i,i) =1 ;
     elseif (y(i) == 0 && x(i) >= 0 && x(i) <= xmax) 
        TMatrix(i,:) = 0;
        TMatrix(i,i) = 1;
        SPM(i,:) = 0;
        LV(i) = 0;
    elseif ( y(i) == xmax && x(i) >= 0 && x(i) <= xmax) 
        LV(i) = 0;
        TMatrix(i,:) = 0;
        TMatrix(i,i) =1 ;
        SPM(i,:) = 0 ; 
    end
end

for j = 1:M+1
    RHS = SPM*Q+LV;
    Q = TMatrix\RHS;
    trisurf(LNODES,x,y,Q(:,:))
shading interp
xlabel('x','fontsize',16) 
xlim([0 xmax])
ylim([0 xmax])
zlim([0 1])
ylabel('y','fontsize',16)
zlabel('q(x,y)','fontsize',16)
title(['Heat diffusion at t= ',num2str(T(j))],'fontsize',20)
pause(1e-10)  
end

max(Q)



