%% This script solves schnakenberg reaction kinetics on a fixed square.
% Wakil Sarfaraz   30/11/15
clear all;
tic
xmax = 1;
tm = 1;
dt = 0.01;
M = tm/dt;

d_1 = 1;
d_2 = 1;
a = 0.1;
b = 0.9;

N = 40; 
X = linspace(0,xmax,N+1);
T = linspace(0, tm,M+1); 
[x, y] = meshgrid(X,X); 
                       
x =x(:);  
y =y(:);  
NNODES = (N+1)^2;
U = zeros(NNODES,1);
V = zeros(NNODES,1);
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

SPU = sparse(NNODES, NNODES);
SPV = sparse(NNODES, NNODES);
SPMU = sparse(NNODES, NNODES);
SPMV = sparse(NNODES, NNODES);
SPA1 = sparse(NNODES,NNODES);
SPA2 =sparse (NNODES,NNODES);
LU = zeros(NNODES,1);         
LV = zeros(NNODES,1);

for n = 1: NTRI
    r1 = [x(LNODES(n,1)) y(LNODES(n,1))];
    r2 = [x(LNODES(n,2)) y(LNODES(n,2))];
    r3 = [x(LNODES(n,3)) y(LNODES(n,3))];
    J = [r2(1)-r1(1) r2(2)-r1(2); r3(1)-r1(1) r3(2)-r1(2)]; 
    
 Astiff =(1/det(J))*eye(3,3);     
 Amass = (1/(120*det(J)))*eye(3,3);   
 A1 = (1/det(J))*eye(3,3);
 A2 = (1/det(J))*eye(3,3);
   
   
   
   
       for i = 1 : 3 
           for j = 1:3  
               SPU(LNODES(n,i),LNODES(n,j)) =SPU(LNODES(n,i),LNODES(n,j))+ Astiff(i,j);
               SPV(LNODES(n,i),LNODES(n,j)) =SPV(LNODES(n,i),LNODES(n,j))+ Astiff(i,j);
               SPMU(LNODES(n,i),LNODES(n,j)) =SPMU(LNODES(n,i),LNODES(n,j))+ Amass(i,j);
               SPMV(LNODES(n,i),LNODES(n,j)) = SPMV(LNODES(n,i),LNODES(n,j))+Amass(i,j);
               SPA1(LNODES(n,i),LNODES(n,j)) =  SPA1(LNODES(n,i),LNODES(n,j))+ A1(i,j);
               SPA2(LNODES(n,i),LNODES(n,j)) =  SPA2(LNODES(n,i),LNODES(n,j))+ A2(i,j);
           end
       end
       
    
       
end


TMatrixU =  ((1/dt)*SPMU-d_1*SPU+SPA1);
TMatrixV = (((dt+1)/dt)*SPMV-d_2*SPV-SPA2);

LU = (1/dt)*SPMU*U+a;
LV = (1/dt)*SPMV*V+b;



for i = 1: NNODES
    if (x(i)==xmax && y(i)>=0 && y(i) <= xmax)  
        LU(i) = 0;
        LV(i) = 0;
        TMatrixU(i,:) = 0;
        TMatrixV(i,:) = 0;
        SPMU(i,:) = 0;
        SPMV(i,:) = 0;
        TMatrixU(i,i) = 1;
        TMatrixV(i,i) = 1;
    elseif (x(i)==0 && y(i)>=0 && y(i) <= xmax)  
        LU(i) = 0;
        LV(i) = 0;
        TMatrixU(i,:) = 0;
        TMatrixV(i,:) = 0;
        SPMU(i,:) = 0;
        SPMV(i,:) = 0;
        TMatrixU(i,i) =1;
        TMatrixV(i,i) =1;
     elseif (y(i) == 0 && x(i) >= 0 && x(i) <= xmax) 
        TMatrixU(i,:) = 0;
        TMatrixV(i,:) = 0;
        TMatrixU(i,i) = 1;
        TMatrixV(i,i) = 1;
        SPMU(i,:) = 0;
        SPMV(i,:) = 0;
        LU(i) = 0;
        LV(i) = 0;
    elseif ( y(i) == xmax && x(i) >= 0 && x(i) <= xmax) 
        TMatrixU(i,:) = 0;
        TMatrixV(i,:) = 0;
        TMatrixU(i,i) = 1;
        TMatrixV(i,i) = 1;
        SPMU(i,:) = 0;
        SPMV(i,:) = 0;
        LU(i) = 0;
        LV(i) = 0;
    end
end

for j = 1:M+1
    RHSU = (1/dt)*SPMU*U+LU+a;
    RHSV = (1/dt)*SPMV*V+LV+b;
    U = TMatrixU\RHSU;
    V = TMatrixV\RHSV;
    trisurf(LNODES,x,y,U(:,:),V(:,:))
%shading interp
xlabel('x','fontsize',16) 
xlim([0 xmax])
ylim([0 xmax])
zlim([0 100])
ylabel('y','fontsize',16)
zlabel('u & v','fontsize',16)
title(['Evolution of reaction kinetics at t= ',num2str(T(j))],'fontsize',16)
pause(1e-10) 

end


 trisurf(LNODES,x,y,U(:,:),V(:,:))
 xlabel('x')
 ylabel('y')
 zlabel('u and v')
 %view(2)
 title ('Schankenberg Kinetics with d1=0.1. d2=0.1')
 %shading interp
 
 
 
 
 
 
 
 