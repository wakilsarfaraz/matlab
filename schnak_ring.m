%% This script solves schnakenberg reaction kinetics on a fixed circle.
% Wakil Sarfaraz   30/11/15
clear all;

addpath distmesh

% xmax = 1;
tm = 1;
dt = 0.01;
M = tm/dt;

du = 0.08;
dv = 0.01;
a = 0.1;
b = 0.9;
gamma = 300000;
T = linspace(0, tm,M+1); 

  fd=inline('-0.1+abs(0.2-sqrt(sum(p.^2,2)))');
  [p,t]=distmesh2d(fd,@huniform,0.01,[-1,-1;1,1],[]);

x = p(:,1);
y = p(:,2);

NNODES = length(x);
NTRI = size(t,1);
LNODES = t;



U = zeros(NNODES,1);
V = zeros(NNODES,1);


for i = 1 : NNODES

    U(i) = U(i)+0.0001*pi^2*sin(100*pi*sin(5*pi*x(i))+10*pi*sin(100*pi*y(i)));%*sin(0.5*pi*x(i));
    %V(i) = V(i)+0.0001*pi^2*exp(-x(i)*y(i))*(sin(15*pi*x(i))+cos(15*pi*y(i)));
    V(i) = V(i)+0.001*sin(50*pi*x(i)+5*pi*y(i));
end
% figure(1)
% title ('Schnakenberg Kinetics with Du=40, Dv=2, gamma=500')
% subplot (2,2,1)
% trisurf(LNODES,x,y,U(:,:))
% xlabel('x')
% ylabel('y')
% view(2)
% legend('Initial u')
% shading interp
% axis equal tight
% subplot(2,2,2)
% trisurf(LNODES,x,y,V(:,:))
% xlabel('x')
% ylabel('y')
% view(2)
% legend('Initial v')
% shading interp
% axis equal tight



SPSM = sparse(NNODES, NNODES);
SPMM = sparse(NNODES, NNODES);
SPC = sparse(NNODES,NNODES);
SPD =sparse (NNODES,NNODES);

UG = zeros(NNODES,1);
VG = zeros(NNODES,1);

RHSU = zeros(NNODES,1);         
RHSV = zeros(NNODES,1);

for n = 1: NTRI
    r1 = [x(LNODES(n,1)) y(LNODES(n,1))];
    r2 = [x(LNODES(n,2)) y(LNODES(n,2))];
    r3 = [x(LNODES(n,3)) y(LNODES(n,3))];
    J = [r2(1)-r1(1) r2(2)-r1(2); r3(1)-r1(1) r3(2)-r1(2)]; 
    
 StiffL = (1/(2*det(J)))* [(r2-r3)*(r2-r3)' (r2-r3)*(r3-r1)' (r2-r3)*(r1-r2)';... 
           (r2-r3)*(r3-r1)' (r3-r1)*(r3-r1)' (r3-r1)*(r1-r2)';...
           (r2-r3)*(r1-r2)' (r3-r1)*(r1-r2)' (r1-r2)*(r1-r2)']; 
       
 MassL = det(J)/120*[1/12 1/24 1/24; 1/24 1/12 1/24; 1/24 1/24 1/12]; 

 
 CL = (det(J))*[(1/15*V(LNODES(n,1))-1/30*V(LNODES(n,2))-1/30*V(LNODES(n,3)))*U(LNODES(n,1))+...
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
 
 
 DL = (det(J))*[1/15*(U(LNODES(n,1))*U(LNODES(n,1))-U(LNODES(n,1))*U(LNODES(n,2))-U(LNODES(n,1))*U(LNODES(n,3)))+1/9*...
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
 
 UL = a*det(J)*[0;1/2;1/2];
 VL = b*det(J)*[0;1/2;1/2];
   
   
   
   
       for i = 1 : 3 
           for j = 1:3  
               SPSM(LNODES(n,i),LNODES(n,j)) =SPSM(LNODES(n,i),LNODES(n,j))+ StiffL(i,j);
               SPMM(LNODES(n,i),LNODES(n,j)) =SPMM(LNODES(n,i),LNODES(n,j))+ MassL(i,j);
               SPC(LNODES(n,i),LNODES(n,j)) =  SPC(LNODES(n,i),LNODES(n,j))+ CL(i,j);
               SPD(LNODES(n,i),LNODES(n,j)) =  SPD(LNODES(n,i),LNODES(n,j))+ DL(i,j);
           end
       end
       
       
       for i = 1 : 3
           UG(LNODES(n,i)) = UG(LNODES(n,i))+UL(i);
           VG(LNODES(n,i)) = VG(LNODES(n,i))+VL(i);
       end
           
       
    
       
end


 TMatrixU =  SPMM-dt*du*SPSM+dt*gamma*SPMM-dt*gamma*SPC;
 TMatrixV =  SPMM-dt*dv*SPSM+gamma*SPD;

for j = 1:M+1


    RHSU = SPMM*U+dt*gamma*UG;
    RHSV = SPMM*V+dt*gamma*VG;
    U = TMatrixU\RHSU;
    V = TMatrixV\RHSV;
    
    figure(1)
    
%subplot(1,2,1)
trisurf(LNODES,x,y,1/max(U)*U(:,:))
shading interp
xlabel('x','fontsize',16) 
% xlim([0 xmax])
% ylim([0 xmax])
view(2)
ylabel('y','fontsize',16)
zlabel('u & v','fontsize',16)
title(['Evolution of u at t= ',num2str(T(j))],'fontsize',8)
axis equal tight
% subplot(1,2,2)
% trisurf(LNODES,x,y,1/max(V)*V(:,:))
% shading interp
% xlabel('x','fontsize',16) 
% % xlim([0 xmax])
% % ylim([0 xmax])
% view(2)
% ylabel('y','fontsize',16)
% zlabel('u & v','fontsize',16)
% title(['Evolution of v at t= ',num2str(T(j))],'fontsize',8)
% axis equal tight
%MV(j)=getframe(gcf);
pause(1e-10) 



end


figure(2)
subplot(1,2,1)
 trisurf(LNODES,x,y,1/max(U)*U(:,:))
 xlabel('x')
 ylabel('y')
 view(2)
 legend('Evolved pattern of u')
 shading interp
 axis equal tight
 subplot(1,2,2)
 trisurf(LNODES,x,y,1/max(V)*V(:,:))
 xlabel('x')
 ylabel('y')
 view(2)
 legend('Evolved pattern of v')
 shading interp
 axis equal tight
%movie2avi(MV,'SpotsToSptripes.avi');