%% This script solves schnakenberg reaction kinetics on a fixed square.
% Wakil Sarfaraz   30/11/15

tic
clear all; close all; clc;
xmax = 1;
tm = 4;
dt = .005;
M = tm/dt;
du = .1;
dv = .1;
a = .1;
b = 1.2;
gam = 1;
eps = .25;
N =64; 
X = linspace(0,xmax,N+1);
T = linspace(0,tm,M+1); 
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
    
        U(i) = a + b +0.05*(cos(6*pi*(x(i)))*cos(6*pi*y(i)));
    V(i) = b/(a+b)^2+0.05*(cos(6*pi*x(i))*cos(6*pi*x(i)));
% 
%     U(i) = a + b +0.001*(cos(4*pi*(x(i)))*cos(4*pi*y(i)));
%     V(i) = b/(a+b)^2+0.001*(cos(4*pi*x(i))*cos(4*pi*x(i)));

% U(i) = 0.919145 + 0.0016*cos(4*pi*(x(i)+y(i)))+0.01*(cos(1*pi*x(i))+...
%     cos(2*pi*x(i))+cos(3*pi*x(i))+cos(4*pi*x(i))+cos(5*pi*x(i))+...
%     cos(6*pi*x(i))+cos(7*pi*x(i))+cos(8*pi*x(i)));
% V(i) = 0.937903 + 0.0016*cos(4*pi*(x(i)+y(i)))+0.01*(cos(1*pi*x(i))+...
%     cos(2*pi*x(i))+cos(3*pi*x(i))+cos(4*pi*x(i))+cos(5*pi*x(i))+...
%     cos(6*pi*x(i))+cos(7*pi*x(i))+cos(8*pi*x(i)));

% U(i) = a + b + 0.0016*cos(6*pi*(x(i)+y(i)))+.001*(cos(1*pi*x(i))+...
%     cos(2*pi*x(i))+cos(3*pi*x(i))+cos(4*pi*x(i))+cos(5*pi*x(i))+...
%     cos(6*pi*x(i))+cos(7*pi*x(i))+cos(8*pi*x(i)));
% V(i) = b/(a+b)^2 + 0.0016*cos(6*pi*(x(i)+y(i)))+.001*(cos(1*pi*x(i))+...
%     cos(2*pi*x(i))+cos(3*pi*x(i))+cos(4*pi*x(i))+cos(5*pi*x(i))+...
%     cos(6*pi*x(i))+cos(7*pi*x(i))+cos(8*pi*x(i)));

end
SPSM = sparse(NNODES, NNODES);
SPMM = sparse(NNODES, NNODES);
SPC = sparse(NNODES,NNODES);
SPD =sparse (NNODES,NNODES);

A = zeros(NNODES,1);


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
 AL = det(J)/6*[1;1;1];
       for i = 1 : 3 
           for j = 1:3  
               SPSM(LNODES(n,i),LNODES(n,j)) =SPSM(LNODES(n,i),LNODES(n,j))+ StiffL(i,j);
               SPMM(LNODES(n,i),LNODES(n,j)) =SPMM(LNODES(n,i),LNODES(n,j))+ MassL(i,j);
               SPC(LNODES(n,i),LNODES(n,j)) =  SPC(LNODES(n,i),LNODES(n,j))+ CL(i,j);
               SPD(LNODES(n,i),LNODES(n,j)) =  SPD(LNODES(n,i),LNODES(n,j))+ DL(i,j);
           end
       end
       for l = 1 : 3
           A(LNODES(n,l)) = A(LNODES(n,l))+AL(l);
       end   
end



% for i = 1: NNODES
%     if (x(i)==0 || x(i)==xmax || y(i)==0 || y(i)==xmax)
%         RHSU(i) = 1;
%         RHSV(i) = 1;
%         TMatrixU(i,:) = 0;
%         TMatrixV(i,:) = 0;
%         TMatrixU(i,i) = 1;
%         TMatrixV(i,i) = 1;
%     end
% 
% end

%Neumann B.C
for i = 1 : NTRI
    if(x(LNODES(i,1))==0 || x(LNODES(i,2))==0 ||  x(LNODES(i,3))==0)
        U(LNODES(i,3))=U(LNODES(i,1));
    elseif(x(LNODES(i,1))==xmax || x(LNODES(i,2))==xmax ||  x(LNODES(i,3))==xmax)
        U(LNODES(i,1))=U(LNODES(i,3));
    elseif(y(LNODES(i,1))==0 || y(LNODES(i,2))==0 ||  y(LNODES(i,3))==0)
        U(LNODES(i,2))=U(LNODES(i,1));
    elseif(y(LNODES(i,1))==xmax || y(LNODES(i,2))==xmax ||  y(LNODES(i,3))==xmax)
        U(LNODES(i,1))=U(LNODES(i,2));
    end
end

 TMatrixU =  SPMM+dt*du*SPSM;
 TMatrixV =  SPMM+dt*dv*SPSM;
MatrixU = zeros(length(T),length(U));
MatrixV = zeros(length(T),length(V));
for j = 1:M+1
if(T(j)>=tm/5-eps && T(j)<=tm/5+eps)
    a = b*T(j);
      b = a*T(j);
 gam =.03*gam*T(j);
%     gam = 4*T(j)*gam*exp(abs((tm/4-T(j))^2));
    RHSU = (SPMM-dt*gam*SPMM+dt*gam*SPC)*U+dt*gam*a*A;
    RHSV = (SPMM-gam*dt*SPD)*V+dt*gam*b*A;
    U = TMatrixU\RHSU;
    V = TMatrixV\RHSV;
    MatrixU(j,:)= U;
    MatrixV(j,:)= V;
end
    
% else if (T(j)>=tm/2-eps && T(j)<=tm/2+eps)
% %     a = a*T(j);
% %     b = b*exp((T(j)-tm/2)^2);
% %      gam = gam*exp(-((T(j)-tm/2-eps)^2)/T(j));
% gam = gam*T(j);
%     RHSU = SPMM*U+dt*gam*a*A;
%     RHSV = SPMM*V+dt*gam*b*A;
%     U = TMatrixU\RHSU;
%     V = TMatrixV\RHSV;
%     MatrixU(j,:)= U;
%     MatrixV(j,:)= V;
% else if (T(j)>=3*tm/4-eps && T(j)<=3*tm/4+eps)
% %     a = a*T(j);
% %     b = b*exp(-(T(j)-3*tm/4)^2);
% %     gam = T(j)*gam/exp(-abs((T(j)-3*tm/4)^2));
% gam = gam/T(j);
%     RHSU = SPMM*U+dt*gam*a*A;
%     RHSV = SPMM*V+dt*gam*b*A;
%     U = TMatrixU\RHSU;
%     V = TMatrixV\RHSV;
%     MatrixU(j,:)= U;
%     MatrixV(j,:)= V;
% %  
%     end
%     end
% end
% a = 0.1;
% b = 0.9;

    RHSU = (SPMM-dt*gam*SPMM+dt*gam*SPC)*U+dt*gam*a*A;
    RHSV = (SPMM-gam*dt*SPD)*V+dt*gam*b*A;
    U = TMatrixU\RHSU;
    V = TMatrixV\RHSV;
    MatrixU(j,:)=U;
    MatrixV(j,:)=V;
% figure(1)
% trisurf(LNODES,x,y,U(:,:))
% % colorbar
% shading interp
% xlabel('x','fontsize',20) 
% xlim([0 xmax])
% ylim([0 xmax])
% view(2)
% ylabel('y','fontsize',20)
% zlabel('u & v','fontsize',16)
% title(['Evolution of u at t= ',num2str(T(j))],'fontsize',18)
% title(['Evolution of u at t = 0.1'],'fontsize',18)
% axis equal tight



pause(1e-10) 

end

% 
% figure(1)
% trisurf(LNODES,x,y,1/abs(max(U)-200)*U)
% title('Solution U at t = 10','fontsize',18)
% colorbar
% shading interp
% xlabel('x','fontsize',20) 
% xlim([0 xmax])
% ylim([0 xmax])
% view(2)
% ylabel('y','fontsize',20)
% axis equal tight
%   set(findobj('type','legend'),'fontsize',18)
% set(findobj('type','axes'),'fontsize',18)
 

%figure(2)
L2U = zeros(length(T),1);
L2V = zeros(length(T),1);

for k = 1: length(T)-1
   if (k>=1 && k<=30)
    L2U(k) = sqrt(sum(abs(MatrixU(k+1,:)-MatrixU(k,:)).^2)/(T(j)-T(j-1)));
    L2V(k) = sqrt(sum(abs(MatrixV(k+1,:)-MatrixV(k,:)).^2)/(T(j)-T(j-1)));
   end
    L2U(k) = sqrt(sum(abs(MatrixU(k+1,:)-MatrixU(k,:)).^2)/(T(j)-T(j-1)));
    L2V(k) = sqrt(sum(abs(MatrixV(k+1,:)-MatrixV(k,:)).^2)/(T(j)-T(j-1)));
end

L2U(length(T))=L2U(length(T)-1);
L2V(length(T))=L2V(length(T)-1);
figure(1)
plot(T,L2U,'LineWidth',2,'color','r')

hold on
plot(T,L2V,'LineWidth',2,'color','b')
set(findobj('type','legend'),'fontsize',20)
set(findobj('type','axes'),'fontsize',20)
legend('log(||U^{m+1}-U^m||/\tau)','log(||V^{m+1}-V^m||/\tau)','Location','NorthEast')
xlabel('Time','fontsize',20)
ylabel('log(||\cdot||_L_2)','fontsize',20)
title('Convergence of solutions','fontsize',16)

toc