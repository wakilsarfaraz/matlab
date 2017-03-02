%% This script solves schnakenberg reaction kinetics on a fixed circle.
% Wakil Sarfaraz   30/11/15
clear all;

%addpath distmesh
tic

% tm = 60;
% dt = 0.1;
tm = 9;
dt = .001;

M = tm/dt;

du = .0001;
dv = .01;
% dv = 5;
a = 0.05;
b = 0.9;
gamma = 1000;
T = linspace(0, tm,M+1); 

xmax = 1;
N = 30;
 fd=inline('sqrt(sum(p.^2,2))-1','p');
 [p,t]=distmesh2d(fd,@huniform,xmax/N,[-1,-1;1,1],[]);

x = p(:,1);
y = p(:,2);

[rr tt]=cart2pol(x,y);

NNODES = length(x);
NTRI = size(t,1);
LNODES = t;



U = zeros(NNODES,1);
V = zeros(NNODES,1);


for i = 1 : NNODES
%     if (x(i)==0 || x(i)==xmax || y(i)==0 || y(i)==xmax)
%         U(i) = 0;

%     else
%     U(i) = a + b + 0.1*exp(-10*(abs(sin((x(i)-0.5)^2+(y(i)-1/3)^2))));
%     V(i) = b/(a+b)^2;% + 0.001*exp(sin(-10*((x(i)-0.5)^2+(y(i)-1/3)^2)));
%   U(i) = a + b + 0.05*(cos(6*pi*x(i))+cos(6*pi*y(i)));
%     V(i) = b/(a+b)^2;%+0.1*(sin((x(i)+y(i)))); 
%       U(i) = (cos(2*pi*(x(i))+y(i)))+2.*sin(2*pi*(x(i)^2+y(i)^2));       %eigenmode 9
% U(i) = a+b+(cos(1.8*pi*x(i))*cos(1.8*pi*y(i)))+sin(pi*(x(i)^2+y(i)^2));
U(i) = a+b+cos(4*pi*(x(i)^2))*cos(4*pi*y(i)^2)+0.3*cos(1.5*pi*(x(i)^2+y(i)^2));
% U(i)=a+b+cos(4*pi*(x(i)^2+y(i)^2))+0.3*cos(2*pi*x(i))*cos(2*pi*y(i));
V(i) = b/(a+b)^2;%+cos(2.5*pi*x(i))+cos(2.5*pi*y(i))+(x(i)^2+y(i)^2);
%     V(i) = b/(a+b)^2;%+0.1*(sin((x(i)+y(i)))); 
%     end
end
% for i = 1 : NNODES
% %     if (x(i)==0 || x(i)==xmax || y(i)==0 || y(i)==xmax)
% %         U(i) = 0;
% %         V(i) = 0;
% %     else
% %     U(i) = a + b + exp(-cos(17*(x(i)))-cos(17*y(i)));
%     U(i) = a + b+exp(-cos(50*pi*(x(i)^2+y(i)^2)));%exp(-0.3*cos(3*pi*abs(y(i)))^2-0.3*cos(3*pi*abs(x(i)))^2);%+cos(2*pi*(y(i)^2))+cos(2*pi*x(i)^2);%+ 0.01*exp((((x(i)-0.5)^2+(y(i)-1/3)^2))); 
%     V(i) = b/(a+b)^2;%+exp(-sin(50*pi*(x(i)^2+y(i)^2)));
% %   U(i) = a + b + 0.1*(sin((xmax*pi*x(i)))+sin((xmax*pi*y(i))));
% %   V(i) = b/(a+b)^2;%+0.001*(sin((x(i)+y(i)))); 
% %     end
% end
ui = [min(U) max(U)]
vi = [min(V) max(V)]
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
SPC  = sparse (NNODES,NNODES);
SPD  = sparse (NNODES,NNODES);

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
    
 
 UL = det(J)*[1/6;1/6;1/6];
 VL = det(J)*[1/6;1/6;1/6];
   
   
   
   
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


 TMatrixU =  SPMM+dt*du*SPSM+dt*gamma*SPMM+dt*gamma*SPC;% TMatrixU =  SPMM-dt*du*SPSM+dt*gamma*SPMM-dt*gamma*SPC;
 TMatrixV =  SPMM+dt*dv*SPSM+gamma*SPD;


% for i = 1 : NNODES
%     if (abs(fd(p(i,:)))<=1e-8 )
%         RHSU(i) = 0;
%         RHSV(i) = 0;
%         TMatrixU(i,:) = 0;
%         TMatrixV(i,:) = 0;
%         TMatrixU(i,i) = 1;
%         TMatrixV(i,i) = 1;
%     end
% end
for i = NNODES
    if (abs(fd(p(i,:)))<=1e-8)
        U(LNODES(i,1))=U(LNODES(i,2));
        U(LNODES(i,3))=U(LNODES(i,1));
        U(LNODES(i,2))=U(LNODES(i,3));
        V(LNODES(i,1))=V(LNODES(i,2));
        V(LNODES(i,3))=V(LNODES(i,1));
        V(LNODES(i,2))=V(LNODES(i,3));
    end
end
% Tdiffu = zeros(1,length(T));
% Tdiffv = zeros(1,length(T));
MatrixU = zeros(length(T),length(U));
MatrixV = zeros(length(T),length(U));
for j = 1:M+1




    RHSU = SPMM*U+dt*gamma*a*UG;
    RHSV = SPMM*V+dt*gamma*b*VG;
    U = TMatrixU\RHSU;
    V = TMatrixV\RHSV;
    
%    Tdiffu(j) = sum((U(j+1)-U(j)).^2);
%     Tdiffv(j) = sum((V(j+1)-V(j)).^2);
%     figure(1)
   MatrixU(j,:) = U;
   MatrixV(j,:) = V;
% subplot(1,2,1)
%trisurf(LNODES,x,y,U(:,:))
trisurf(t,x,y,1/max(U)*U(:,:))
 colorbar
shading interp
xlabel('x','fontsize',16) 
view(2)
ylabel('y','fontsize',16)
zlabel('u & v','fontsize',16)
% title(['Evolution of u at t= ',num2str(T(j))],'fontsize',8)
% title('Eigenmode for k=12')
title(['Evolution of u at t= 6'])
axis equal tight
% subplot(1,2,2)
% trisurf(LNODES,x,y,U(:,:))
% %colorbar
% shading interp
% xlabel('x','fontsize',16) 
% view(2)
% ylabel('y','fontsize',16)
% zlabel('u & v','fontsize',16)
% title(['Pattern formed by u] at t= ',num2str(T(j))],'fontsize',8)
% axis equal tight
% %MV(j)=getframe(gcf);
 pause(1e-10) 



end
toc
size(LNODES)
size(x)
% figure(2)
% subplot(1,2,1)
% plot(T,Tdiffu)
% xlabel('Time')
% ylabel('L2 Difference')
% title('U')
% subplot(1,2,2)
% plot(T,Tdiffv)
% title('V')
% xlabel('Time')
% ylabel('L2 Difference')
%figure(2)
% subplot(2,2,3)
%  trisurf(LNODES,x,y,1/max(U)*U(:,:))
%  xlabel('x')
%  ylabel('y')
%  view(2)
%  legend('Evolved pattern of u')
%  shading interp
%  axis equal tight
%  subplot(2,2,4)
%  trisurf(LNODES,x,y,abs(V(:,:)))
%  colorbar
%  xlabel('x')
%  ylabel('y')
%  view(2)
%  legend('Evolved pattern of v')
%  shading interp
%  axis equal tight
%movie2avi(MV,'Disc_pattern.avi');
% figure(1)
% L2U = zeros(length(T),1);
% L2V = zeros(length(T),1);
% 
% for k = 1: length(T)-1
%    
%     L2U(k) = sum(abs(MatrixU(k+1,:)-MatrixU(k,:)).^2);
%     L2V(k) = sum(abs(MatrixV(k+1,:)-MatrixV(k,:)).^2);
% end


% plot(T,L2U,'LineWidth',3,'color','r')
% 
% hold on
% plot(T,L2V,'LineWidth',3,'color','b')
%   set(findobj('type','legend'),'fontsize',18)
% set(findobj('type','axes'),'fontsize',18)
% legend('||U^{m+1}-U^m||_L_2','||V^{m+1}-V^m||_{L_2}')
% xlabel('Time','fontsize',18)
% ylabel('||\cdot||_L_2','fontsize',18)
% title('Convergence of solutions','fontsize',16)
% uf = [min(U) max(U)]
% vf = [min(V) max(V)]