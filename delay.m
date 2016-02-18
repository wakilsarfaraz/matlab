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

% Dh = 0.01;
% Dp = 0.025;

Dh = 0;
Dp = 0;

a1 = 1.2;
a2 = 0.3;
b = 0.2;
c1 = 0.9;
c2 = 0.42;
k1 = 1;
k2 = 2.1;

tau = 0.5; % Delay Parameter

N = 40;
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



for i = 1 : NNODES
    if (x(i)==0 || x(i)==xmax || y(i)==0 || y(i)==xmax)
        H(i)=0;
        P(i)=0;
    else
        H(i) = H(i)+0.0001*pi^2*sin(10*xmax*pi*sin(5*pi*x(i))+10*xmax*pi*sin(5*pi*y(i)));
        P(i) = P(i)+0.0001*pi^2*sin(15*pi*x(i))*cos(15*pi*y(i));
    end
end
figure(1)
subplot(1,2,1)
trisurf(LNODES,x,y,H(:,:))
shading interp
view(2)
axis equal tight
subplot(1,2,2)
trisurf(LNODES,x,y,P(:,:))
shading interp
view(2)
axis equal tight

MG = sparse(NNODES,NNODES);
SG = sparse(NNODES,NNODES);
CG = sparse(NNODES,NNODES);
EG = sparse(NNODES,NNODES);
FG = sparse(NNODES,NNODES);

RH = zeros(NNODES,1);         
RP = zeros(NNODES,1);

for n = 1: NTRI
    r1 = [x(LNODES(n,1)) y(LNODES(n,1))];
    r2 = [x(LNODES(n,2)) y(LNODES(n,2))];
    r3 = [x(LNODES(n,3)) y(LNODES(n,3))];
    J = [r2(1)-r1(1) r2(2)-r1(2); r3(1)-r1(1) r3(2)-r1(2)]; 
    
 SL = (1/(32*det(J)))* [(r2-r3)*(r2-r3)' (r2-r3)*(r3-r1)' (r2-r3)*(r1-r2)';... 
           (r2-r3)*(r3-r1)' (r3-r1)*(r3-r1)' (r3-r1)*(r1-r2)';...
           (r2-r3)*(r1-r2)' (r3-r1)*(r1-r2)' (r1-r2)*(r1-r2)'];      
 ML = det(J)*[1/12 1/24 1/24; 1/24 1/12 1/24; 1/24 1/24 1/12]; 
 CL = (det(J))*k1*[(1/15*V(LNODES(n,1))-1/30*V(LNODES(n,2))-1/30*V(LNODES(n,3)))*U(LNODES(n,1))+...
     (1/18*V(LNODES(n,3))+11/180*V(LNODES(n,2))-1/30*V(LNODES(n,1)))*U(LNODES(n,2))+(11/180*V(LNODES(n,3))+...
     1/18*V(LNODES(n,2))-1/30*V(LNODES(n,1)))*U(LNODES(n,3)) (1/18*V(LNODES(n,3))+11/180*V(LNODES(n,2))-...
     1/30*V(LNODES(n,1)))*U(LNODES(n,1))+(11/180*V(LNODES(n,1))-3/40*V(LNODES(n,2))-5/72*V(LNODES(n,3)))*...
     U(LNODES(n,2))+(1/18*V(LNODES(n,1))-5/72*V(LNODES(n,2))-5/72*V(LNODES(n,3)))*U(LNODES(n,3)) (11/180*...
     V(LNODES(n,3))+1/18*V(LNODES(n,2))-1/30*V(LNODES(n,1)))*U(LNODES(n,1))+(1/18*V(LNODES(n,1))-5/72*V(LNODES(n,2))-...
     5/72*V(LNODES(n,3)))*U(LNODES(n,2))+(11/180*V(LNODES(n,1))-5/72*V(LNODES(n,2))-3/40*V(LNODES(n,3)))*U(LNODES(n,3));...
     (1/18*V(LNODES(n,3))+11/180*V(LNODES(n,2))-1/30*V(LNODES(n,1)))*U(LNODES(n,1))+(11/180*V(LNODES(n,1))-...
     3/40*V(LNODES(n,2))-5/72*V(LNODES(n,3)))*U(LNODES(n,2))+(1/18*V(LNODES(n,1))-5/72*V(LNODES(n,2))-...
     5/72*V(LNODES(n,3)))*U(LNODES(n,3)) (11/180*V(LNODES(n,1))-3/40*V(LNODES(n,2))-5/72*V(LNODES(n,3)))*U(LNODES(n,1))+...
     (1/8*V(LNODES(n,3))+1/5*V(LNODES(n,2))-3/40*V(LNODES(n,1)))*U(LNODES(n,2))+(1/9*V(LNODES(n,3))+1/8*V(LNODES(n,2))-...
     5/72*V(LNODES(n,1)))*U(LNODES(n,3)) (1/18*V(LNODES(n,1))-5/72*V(LNODES(n,2))-5/72*V(LNODES(n,3)))*U(LNODES(n,1))+...
     (1/8*V(LNODES(n,3))+1/9*V(LNODES(n,2))-5/72*V(LNODES(n,1)))*U(LNODES(n,2))+(1/5*V(LNODES(n,3))+1/8*V(LNODES(n,2))-...
     3/40*V(LNODES(n,1)))*U(LNODES(n,3)); (11/180*V(LNODES(n,3))+1/18*V(LNODES(n,2))-1/30*V(LNODES(n,1)))*U(LNODES(n,1))+...
     (1/18*V(LNODES(n,1))-5/72*V(LNODES(n,2))-5/72*V(LNODES(n,3)))*U(LNODES(n,2))+(11/180*V(LNODES(n,1))-5/72*V(LNODES(n,2))-...
     3/40*V(LNODES(n,3)))*U(LNODES(n,3)) (1/18*V(LNODES(n,1))-5/72*V(LNODES(n,2))-5/72*V(LNODES(n,3)))*U(LNODES(n,1))+...
     (1/8*V(LNODES(n,3))+1/9*V(LNODES(n,2))-5/72*V(LNODES(n,1)))*U(LNODES(n,2))+(1/5*V(LNODES(n,3))+1/8*V(LNODES(n,2))-...
     3/40*V(LNODES(n,1)))*U(LNODES(n,3)) (11/180*V(LNODES(n,1))-5/72*V(LNODES(n,2))-3/40*V(LNODES(n,3)))*U(LNODES(n,1))+...
     (1/8*V(LNODES(n,3))+1/9*V(LNODES(n,2))-5/72*V(LNODES(n,1)))*U(LNODES(n,2))+(1/5*V(LNODES(n,3))+1/8*V(LNODES(n,2))-...
     3/40*V(LNODES(n,1)))*U(LNODES(n,3))];
 EL = (det(J))*k2*[1/15*(U(LNODES(n,1))*U(LNODES(n,1))-U(LNODES(n,1))*U(LNODES(n,2))-U(LNODES(n,1))*U(LNODES(n,3)))+1/9*...
     U(LNODES(n,2))*U(LNODES(n,3))+11/180*(U(LNODES(n,2))*U(LNODES(n,2))+U(LNODES(n,3))*U(LNODES(n,3))) 1/9*U(LNODES(n,1))*...
     U(LNODES(n,3))+11/90*U(LNODES(n,1))*U(LNODES(n,2))-1/30*U(LNODES(n,1))*U(LNODES(n,1))-5/36*U(LNODES(n,2))*U(LNODES(n,3))-...
     3/40*U(LNODES(n,2))*U(LNODES(n,2))-5/72*U(LNODES(n,3))*U(LNODES(n,3)) 11/90*U(LNODES(n,1))*U(LNODES(n,3))+1/9*U(LNODES(n,1))*...
     U(LNODES(n,2))-1/30*U(LNODES(n,1))*U(LNODES(n,1))-5/36*U(LNODES(n,2))*U(LNODES(n,3))-5/72*U(LNODES(n,2))*U(LNODES(n,2))-...
     3/40*U(LNODES(n,3))*U(LNODES(n,3));1/9*U(LNODES(n,1))*U(LNODES(n,3))+11/90*U(LNODES(n,1))*U(LNODES(n,2))-1/30*...
     U(LNODES(n,1))*U(LNODES(n,1))-5/36*U(LNODES(n,2))*U(LNODES(n,3))-3/40*U(LNODES(n,2))*U(LNODES(n,2))-5/72*...
     U(LNODES(n,3))*U(LNODES(n,3)) 11/180*U(LNODES(n,1))*U(LNODES(n,1))-3/20*U(LNODES(n,2))*U(LNODES(n,1))-5/36*U(LNODES(n,3))*...
     U(LNODES(n,1))+1/4*U(LNODES(n,3))*U(LNODES(n,2))+1/5*U(LNODES(n,2))*U(LNODES(n,2))+1/9*U(LNODES(n,3))*U(LNODES(n,3)) 1/18*...
     U(LNODES(n,1))*U(LNODES(n,1))-5/36*U(LNODES(n,2))*U(LNODES(n,1))-5/36*U(LNODES(n,3))*U(LNODES(n,1))+2/9*U(LNODES(n,3))*...
     U(LNODES(n,2))+1/8*U(LNODES(n,2))*U(LNODES(n,2))+1/9*U(LNODES(n,3))*U(LNODES(n,3));11/90*U(LNODES(n,1))*U(LNODES(n,3))+1/9*...
     U(LNODES(n,1))*U(LNODES(n,2))-1/30*U(LNODES(n,1))*U(LNODES(n,1))-5/36*U(LNODES(n,2))*U(LNODES(n,3))-5/72*U(LNODES(n,2))*...
     U(LNODES(n,2))-3/40*U(LNODES(n,3))*U(LNODES(n,3)) 1/18*U(LNODES(n,1))*U(LNODES(n,1))-5/36*U(LNODES(n,2))*U(LNODES(n,1))-...
     5/36*U(LNODES(n,3))*U(LNODES(n,1))+2/9*U(LNODES(n,3))*U(LNODES(n,2))+1/8*U(LNODES(n,2))*U(LNODES(n,2))+1/9*U(LNODES(n,3))*...
     U(LNODES(n,3)) 11/180*U(LNODES(n,1))*U(LNODES(n,1))-5/36*U(LNODES(n,1))*U(LNODES(n,2))-3/20*U(LNODES(n,1))*U(LNODES(n,3))+...
     1/9*U(LNODES(n,2))*U(LNODES(n,2))+1/4*U(LNODES(n,2))*U(LNODES(n,2))+1/5*U(LNODES(n,3))*U(LNODES(n,3))];
 FL = (det(J))*[(1/15*V(LNODES(n,1))-1/30*V(LNODES(n,2))-1/30*V(LNODES(n,3)))*U(LNODES(n,1))+...
     (1/18*V(LNODES(n,3))+11/180*V(LNODES(n,2))-1/30*V(LNODES(n,1)))*U(LNODES(n,2))+(11/180*V(LNODES(n,3))+...
     1/18*V(LNODES(n,2))-1/30*V(LNODES(n,1)))*U(LNODES(n,3)) (1/18*V(LNODES(n,3))+11/180*V(LNODES(n,2))-...
     1/30*V(LNODES(n,1)))*U(LNODES(n,1))+(11/180*V(LNODES(n,1))-3/40*V(LNODES(n,2))-5/72*V(LNODES(n,3)))*...
     U(LNODES(n,2))+(1/18*V(LNODES(n,1))-5/72*V(LNODES(n,2))-5/72*V(LNODES(n,3)))*U(LNODES(n,3)) (11/180*...
     V(LNODES(n,3))+1/18*V(LNODES(n,2))-1/30*V(LNODES(n,1)))*U(LNODES(n,1))+(1/18*V(LNODES(n,1))-5/72*V(LNODES(n,2))-...
     5/72*V(LNODES(n,3)))*U(LNODES(n,2))+(11/180*V(LNODES(n,1))-5/72*V(LNODES(n,2))-3/40*V(LNODES(n,3)))*U(LNODES(n,3));...
     (1/18*V(LNODES(n,3))+11/180*V(LNODES(n,2))-1/30*V(LNODES(n,1)))*U(LNODES(n,1))+(11/180*V(LNODES(n,1))-...
     3/40*V(LNODES(n,2))-5/72*V(LNODES(n,3)))*U(LNODES(n,2))+(1/18*V(LNODES(n,1))-5/72*V(LNODES(n,2))-...
     5/72*V(LNODES(n,3)))*U(LNODES(n,3)) (11/180*V(LNODES(n,1))-3/40*V(LNODES(n,2))-5/72*V(LNODES(n,3)))*U(LNODES(n,1))+...
     (1/8*V(LNODES(n,3))+1/5*V(LNODES(n,2))-3/40*V(LNODES(n,1)))*U(LNODES(n,2))+(1/9*V(LNODES(n,3))+1/8*V(LNODES(n,2))-...
     5/72*V(LNODES(n,1)))*U(LNODES(n,3)) (1/18*V(LNODES(n,1))-5/72*V(LNODES(n,2))-5/72*V(LNODES(n,3)))*U(LNODES(n,1))+...
     (1/8*V(LNODES(n,3))+1/9*V(LNODES(n,2))-5/72*V(LNODES(n,1)))*U(LNODES(n,2))+(1/5*V(LNODES(n,3))+1/8*V(LNODES(n,2))-...
     3/40*V(LNODES(n,1)))*U(LNODES(n,3)); (11/180*V(LNODES(n,3))+1/18*V(LNODES(n,2))-1/30*V(LNODES(n,1)))*U(LNODES(n,1))+...
     (1/18*V(LNODES(n,1))-5/72*V(LNODES(n,2))-5/72*V(LNODES(n,3)))*U(LNODES(n,2))+(11/180*V(LNODES(n,1))-5/72*V(LNODES(n,2))-...
     3/40*V(LNODES(n,3)))*U(LNODES(n,3)) (1/18*V(LNODES(n,1))-5/72*V(LNODES(n,2))-5/72*V(LNODES(n,3)))*U(LNODES(n,1))+...
     (1/8*V(LNODES(n,3))+1/9*V(LNODES(n,2))-5/72*V(LNODES(n,1)))*U(LNODES(n,2))+(1/5*V(LNODES(n,3))+1/8*V(LNODES(n,2))-...
     3/40*V(LNODES(n,1)))*U(LNODES(n,3)) (11/180*V(LNODES(n,1))-5/72*V(LNODES(n,2))-3/40*V(LNODES(n,3)))*U(LNODES(n,1))+...
     (1/8*V(LNODES(n,3))+1/9*V(LNODES(n,2))-5/72*V(LNODES(n,1)))*U(LNODES(n,2))+(1/5*V(LNODES(n,3))+1/8*V(LNODES(n,2))-...
     3/40*V(LNODES(n,1)))*U(LNODES(n,3))];
       for i = 1 : 3 
           for j = 1:3  
               MG(LNODES(n,i),LNODES(n,j))  =  MG(LNODES(n,i),LNODES(n,j))+ ML(i,j);
               SG(LNODES(n,i),LNODES(n,j))  =  SG(LNODES(n,i),LNODES(n,j))+ SL(i,j);
               CG(LNODES(n,i),LNODES(n,j))  =  CG(LNODES(n,i),LNODES(n,j))+ CL(i,j);
               EG(LNODES(n,i),LNODES(n,j))  =  EG(LNODES(n,i),LNODES(n,j))+ EL(i,j);
               FG(LNODES(n,i),LNODES(n,j))  =  FG(LNODES(n,i),LNODES(n,j))+ FL(i,j);
           end
       end   
end

TH = (1/dt-a1)*MG-Dh*SG+b*CG;
TP = (1/dt-a2)*MG-Dp*SG+c2*FG;

for i = 1: NNODES
    if (x(i)==0 || x(i)==xmax || y(i)==0 || y(i)==xmax)
        RH(i) = 0;
        RP(i) = 0;
        TH(i,:) = 0;
        TP(i,:) = 0;
        TH(i,i) = 1;
        TP(i,i) = 1;
    end
end

for j = 1 : M+1
    RH = 1/dt*MG*H-c1*EG*P;
    
    H = TH\RH;
    P = TP\RP;
    
    figure(2)
    subplot(1,2,1)
    trisurf(LNODES,x,y,H(:,:));
    shading interp
    xlabel('x','fontsize',16) 
    xlim([0 xmax])
    ylim([0 xmax])
    view(2)
    ylabel('y','fontsize',16)
    zlabel('H','fontsize',16)
    title(['Evolution of H at t= ',num2str(T(j))],'fontsize',8)
    axis equal tight
    subplot(1,2,2)
    trisurf(LNODES,x,y,P(:,:));
    shading interp
    xlabel('x','fontsize',16) 
    view(2)
    ylabel('y','fontsize',16)
    zlabel('P','fontsize',16)
    title(['Evolution of P at t= ',num2str(T(j))],'fontsize',8)
    axis equal tight
    pause(1e-5)
end
    










































































