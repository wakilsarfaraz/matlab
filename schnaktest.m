%% This script solves schnakenberg reaction kinetics on a fixed square.
% Wakil Sarfaraz   30/11/15
clear all;
tic
xmax = 1;
tm = 1;
dt = 0.01;
M = tm/dt;

du = 20;
dv = 15;
a = 0.1;
b = 0.9;
gam = 86;

N = 50; 
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

for i = 1 : NNODES
    if (x(i)==0 || x(i)==xmax || y(i)==0 || y(i)==xmax)
        U(i) = 0;
        V(i) = 0;
    else
    U(i) = a + b + 0.001*exp(-100*(abs(sin((x(i)-0.5)^2+(y(i)-1/3)^2))));
    V(i) = a + b + 0.001*abs(sin(-100*((x(i)-0.5)^2+(y(i)-1/3)^2)));
  
    %V(i) = b/(a+b)^2; 
    end
end
ui = [min(U) max(U)]
vi = [min(V) max(V)]
% figure(1)
% title ('Schnakenberg Kinetics with Du=40, Dv=2, gamma=500')
% subplot (2,2,1)
% trisurf(LNODES,x,y,U(:,:))
% xlabel('x')
% ylabel('y')
% view(2)
% %title('Initial pattern of u')
% legend('Initial u')
% shading interp
% axis equal tight
% subplot(2,2,2)
% trisurf(LNODES,x,y,V(:,:))
% xlabel('x')
% ylabel('y')
% view(2)
% %title('Initial pattern of v')
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
       
       ksi = [0 1/2]';
       eta = [0 1/2]';
       
       xx = [(1-ksi(1)-eta(1))*r1(1)+ksi(1)*r2(1)+eta(1)*r3(1) (1-ksi(2)-eta(2))*r1(1)+ksi(2)*r2(1)+eta(2)*r3(1)]';
       yy = [(1-ksi(1)-eta(1))*r1(2)+ksi(1)*r2(2)+eta(1)*r3(2) (1-ksi(2)-eta(2))*r1(2)+ksi(2)*r2(2)+eta(2)*r3(2)]';
       
 MassL = det(J)*[1/12 1/24 1/24; 1/24 1/12 1/24; 1/24 1/24 1/12]; 
 

 CL = det(J)/360*[12*U(LNODES(n,1))*V(LNODES(n,1))+3*(U(LNODES(n,1))*V(LNODES(n,2))+U(LNODES(n,2))*V(LNODES(n,1))+...
     U(LNODES(n,1))*V(LNODES(n,3))+U(LNODES(n,3))*V(LNODES(n,1)))+2*(U(LNODES(n,2))*V(LNODES(n,2))+U(LNODES(n,3))*V(LNODES(n,3)))+...
     U(LNODES(n,2))*V(LNODES(n,3))+U(LNODES(n,3))*V(LNODES(n,2)) 3*(U(LNODES(n,1))*V(LNODES(n,1))+U(LNODES(n,2))*V(LNODES(n,2)))+...
     2*(U(LNODES(n,1))*V(LNODES(n,2))+U(LNODES(n,2))*V(LNODES(n,1)))+V(LNODES(n,3))*(U(LNODES(n,1))+U(LNODES(n,2))+U(LNODES(n,3)))+...
     U(LNODES(n,3))*(V(LNODES(n,1))+V(LNODES(n,2))) 3*(U(LNODES(n,1))*V(LNODES(n,1))+U(LNODES(n,3))*V(LNODES(n,3)))+2*...
     (U(LNODES(n,1))*V(LNODES(n,3))+U(LNODES(n,3))*V(LNODES(n,1)))+V(LNODES(n,2))*(U(LNODES(n,1))+U(LNODES(n,2))+U(LNODES(n,3)))+...
     U(LNODES(n,2))*(V(LNODES(n,1))+V(LNODES(n,3))); 3*(U(LNODES(n,1))*V(LNODES(n,1))+U(LNODES(n,2))*V(LNODES(n,2)))+...
     2*(U(LNODES(n,1))*V(LNODES(n,2))+U(LNODES(n,2))*V(LNODES(n,1)))+V(LNODES(n,3))*(U(LNODES(n,1))+U(LNODES(n,2))+U(LNODES(n,3)))+...
     U(LNODES(n,3))*(V(LNODES(n,1))+V(LNODES(n,2))) 12*U(LNODES(n,2))*V(LNODES(n,2))+3*(U(LNODES(n,1))*V(LNODES(n,2))+U(LNODES(n,2))*V(LNODES(n,1))+...
     U(LNODES(n,2))*V(LNODES(n,3))+U(LNODES(n,3))*V(LNODES(n,2)))+2*(U(LNODES(n,1))*V(LNODES(n,1))+U(LNODES(n,3))*V(LNODES(n,3)))+...
     U(LNODES(n,1))*V(LNODES(n,3))+U(LNODES(n,3))*V(LNODES(n,1)) 3*(U(LNODES(n,2))*V(LNODES(n,2))+U(LNODES(n,3))*V(LNODES(n,3)))+...
     2*(U(LNODES(n,2))*V(LNODES(n,3))+U(LNODES(n,3))*V(LNODES(n,2)))+V(LNODES(n,1))*(U(LNODES(n,1))+U(LNODES(n,2))+U(LNODES(n,3)))+...
     U(LNODES(n,1))*(V(LNODES(n,2))+V(LNODES(n,3))); 3*(U(LNODES(n,1))*V(LNODES(n,1))+U(LNODES(n,3))*V(LNODES(n,3)))+2*...
     (U(LNODES(n,1))*V(LNODES(n,3))+U(LNODES(n,3))*V(LNODES(n,1)))+V(LNODES(n,2))*(U(LNODES(n,1))+U(LNODES(n,2))+U(LNODES(n,3)))+...
     U(LNODES(n,2))*(V(LNODES(n,1))+V(LNODES(n,3))) 3*(U(LNODES(n,2))*V(LNODES(n,2))+U(LNODES(n,3))*V(LNODES(n,3)))+...
     2*(U(LNODES(n,2))*V(LNODES(n,3))+U(LNODES(n,3))*V(LNODES(n,2)))+V(LNODES(n,1))*(U(LNODES(n,1))+U(LNODES(n,2))+U(LNODES(n,3)))+...
     U(LNODES(n,1))*(V(LNODES(n,2))+V(LNODES(n,3))) 12*U(LNODES(n,3))*V(LNODES(n,3))+3*(U(LNODES(n,1))*V(LNODES(n,3))+...
     U(LNODES(n,3))*V(LNODES(n,1))+U(LNODES(n,2))*V(LNODES(n,3))+U(LNODES(n,3))*V(LNODES(n,2)))+2*(U(LNODES(n,1))*V(LNODES(n,1))+...
     U(LNODES(n,2))*V(LNODES(n,2)))+U(LNODES(n,1))*V(LNODES(n,2))+U(LNODES(n,2))*V(LNODES(n,1))];
 

DL = det(J)/360*[12*U(LNODES(n,1))^2+6*(U(LNODES(n,1))*U(LNODES(n,2))+U(LNODES(n,1))*U(LNODES(n,3)))+2*(U(LNODES(n,2))^2+...
    U(LNODES(n,3))^2+U(LNODES(n,2))*U(LNODES(n,3))) 4*U(LNODES(n,1))*U(LNODES(n,2))+3*(U(LNODES(n,1))^2+U(LNODES(n,2))^2)+...
    2*(U(LNODES(n,1))*U(LNODES(n,3))+U(LNODES(n,2))*U(LNODES(n,3)))+U(LNODES(n,3))^2 4*U(LNODES(n,1))*U(LNODES(n,3))+...
    3*(U(LNODES(n,1))^2+U(LNODES(n,3))^2)+2*(U(LNODES(n,1))*U(LNODES(n,2))+U(LNODES(n,2))*U(LNODES(n,3)))+...
    U(LNODES(n,2))^2; 4*U(LNODES(n,1))*U(LNODES(n,2))+3*(U(LNODES(n,1))^2+U(LNODES(n,2))^2)+...
    2*(U(LNODES(n,1))*U(LNODES(n,3))+U(LNODES(n,2))*U(LNODES(n,3)))+U(LNODES(n,3))^2 12*U(LNODES(n,2))^2+...
    6*(U(LNODES(n,1))*U(LNODES(n,2))+U(LNODES(n,2))*U(LNODES(n,3)))+2*(U(LNODES(n,1))^2+U(LNODES(n,3))^2+...
    U(LNODES(n,1))*U(LNODES(n,3))) 4*U(LNODES(n,2))*U(LNODES(n,3))+3*(U(LNODES(n,2))^2+U(LNODES(n,3))^2)+...
    2*(U(LNODES(n,1))*U(LNODES(n,2))+U(LNODES(n,1))*U(LNODES(n,3)))+U(LNODES(n,1))^2; 4*U(LNODES(n,1))*U(LNODES(n,3))+...
    3*(U(LNODES(n,1))^2+U(LNODES(n,3))^2)+2*(U(LNODES(n,1))*U(LNODES(n,2))+U(LNODES(n,2))*U(LNODES(n,3)))+...
    U(LNODES(n,2))^2 4*U(LNODES(n,2))*U(LNODES(n,3))+3*(U(LNODES(n,2))^2+U(LNODES(n,3))^2)+...
    2*(U(LNODES(n,1))*U(LNODES(n,2))+U(LNODES(n,1))*U(LNODES(n,3)))+U(LNODES(n,1))^2 12*U(LNODES(n,3))^2+...
    6*(U(LNODES(n,1))*U(LNODES(n,3))+U(LNODES(n,2))*U(LNODES(n,3)))+2*(U(LNODES(n,1))^2+U(LNODES(n,2))^2+...
    U(LNODES(n,1))*U(LNODES(n,2)))];
    
 
 UL = a*det(J)*[1/6;1/6;1/6];
 VL = b*det(J)*[1/6;1/6;1/6];
   
   
   
   
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


 TMatrixU =  SPMM-dt*du*SPSM+dt*gam*SPMM-dt*gam*SPC;
 TMatrixV =  SPMM-dt*dv*SPSM+gam*SPD;

for i = 1: NNODES
    if (x(i)==0 || x(i)==xmax || y(i)==0 || y(i)==xmax)
        RHSU(i) = 0;
        RHSV(i) = 0;
        TMatrixU(i,:) = 0;
        TMatrixV(i,:) = 0;
        TMatrixU(i,i) = 1;
        TMatrixV(i,i) = 1;
    end

end

for j = 1:M+1

    RHSU = SPMM*U+dt*gam*UG;
    RHSV = SPMM*V+dt*gam*VG;
    U = TMatrixU\RHSU;
    V = TMatrixV\RHSV;
    
   % figure(1)
    
%subplot(1,2,1)
% trisurf(LNODES,x,y,U(:,:))
% colorbar
% shading interp
% xlabel('x','fontsize',16) 
% xlim([0 xmax])
% ylim([0 xmax])
% view(2)
% ylabel('y','fontsize',16)
% zlabel('u & v','fontsize',16)
% title(['Evolution of u at t= ',num2str(T(j))],'fontsize',8)
% axis equal tight
% subplot(1,2,2)
trisurf(LNODES,x,y,abs(V(:,:)))
colorbar
shading interp
xlabel('x','fontsize',16) 
xlim([0 xmax])
ylim([0 xmax])
view(2)
ylabel('y','fontsize',16)
zlabel('u & v','fontsize',16)
title(['Evolution of v at t= ',num2str(T(j))],'fontsize',8)
axis equal tight
%MV(j)=getframe(gcf);
pause(1e-10) 



end


%figure(2)
% subplot(2,2,3)
%  %trisurf(LNODES,x,y,1/max(U)*U(:,:),1/max(V)*V(:,:))
%  trisurf(LNODES,x,y,U(:,:))
% %  %trisurf(LNODES,x,y,1/max(V)*V(:,:))
%  xlim([0 xmax])
%  ylim([0 xmax])
%  xlabel('x')
%  ylabel('y')
%  %zlabel('u and v')
%  view(2)
%  %title ('Pattern formed by u')
%  legend('Evolved pattern of u')
%  shading interp
%  axis equal tight
%  hold on
%  subplot(2,2,4)
%  %trisurf(LNODES,x,y,1/max(U)*U(:,:),1/max(V)*V(:,:))
%  %trisurf(LNODES,x,y,1/max(U)*U(:,:))
%  trisurf(LNODES,x,y,V(:,:))
%  colorbar 
%  xlim([0 xmax])
%  ylim([0 xmax])
%  xlabel('x')
%  ylabel('y')
%  view(2)
 %title ('Pattern formed by v')
%  legend('Evolved pattern of v')
%  shading interp
%  axis equal tight
%movie2avi(MV,'SpotsToSptripes.avi');
uf = [min(U) max(U)]
vf = [min(V) max(V)]

hold off
 
 
 
 
 
 
 