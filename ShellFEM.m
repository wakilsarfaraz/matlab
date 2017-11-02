%% Wakil Sarfaraz   30/11/15
clear all;

%addpath distmesh
%format long
format short
N = 40;
xmax = 1;

tm = 1;
dt = 0.01;
M = tm/dt;

du = 0.001;
dv = 1;
a = 0.1;
b = 0.9;
epsilon = 0.1;
% gam = 39.596;
gam =.01;
T = linspace(0, tm,M+1); 
  fd=inline('-0.25+abs(0.75-sqrt(sum(p.^2,2)))');
  [p,t]=distmesh2d(fd,@huniform,xmax/N,[-1,-1;1,1],[]);

x = p(:,1);
y = p(:,2);

NNODES = length(x);
NTRI = size(t,1);
LNODES = t;



U = zeros(NNODES,1);
V = zeros(NNODES,1);

for i = 1 : NNODES
      U(i) = a + b + 0.4*(cos(7*pi*x(i))*cos(7*pi*y(i)))+.2*cos(8*pi*(x(i)^2+y(i)^2));
    V(i) = b/(a+b)^2+0.3*(cos(5*pi*(x(i)+y(i)))); 
end



SPSM = sparse(NNODES, NNODES);
SPMM = sparse(NNODES, NNODES);
SPC = sparse(NNODES,NNODES);
SPD =sparse (NNODES,NNODES);

UG = zeros(NNODES,1);
VG = zeros(NNODES,1);

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
       for l = 1 : 3
           A(LNODES(n,l)) = A(LNODES(n,l))+AL(l);
       end    
       
    
       
end
 TMatrixU =  SPMM+dt*du*SPSM+dt*gam*SPMM+dt*gam*SPC;% (original) TMatrixU =  SPMM-dt*du*SPSM+dt*gam*SPMM-dt*gam*SPC;
 TMatrixV =  SPMM+dt*dv*SPSM+gam*SPD;
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



% Tdiffu = zeros(1,length(T));
% Tdiffv = zeros(1,length(T));

MatrixU = zeros(length(T),length(U));
MatrixV = zeros(length(T),length(V));


for j = 1:M+1
    RHSU = (SPMM-dt*gam*SPMM+dt*gam*SPC)*U+dt*gam*a*A;
    RHSV = (SPMM-gam*dt*SPD)*V+dt*gam*b*A;
    U = TMatrixU\RHSU;
    V = TMatrixV\RHSV;
    

% if(T(j)>=tm/30-epsilon && T(j)<=tm/30+epsilon)
%     a = (1/(T(j)-tm/7+eps))*a;
%       b = (1/(T(j)-tm/7)+eps)*b;
%     a = .1*a;
% b = .5*b;
%     RHSU = (SPMM-dt*gam*SPMM+dt*gam*SPC)*U+dt*gam*a*A;
%     RHSV = (SPMM-gam*dt*SPD)*V+dt*gam*b*A;
%     U =TMatrixU\RHSU;
%     V = TMatrixV\RHSV;
%     MatrixU(j,:)=U;
%     MatrixV(j,:)= V;
%     V = 1/min(V)*V;
%     U = 1/min(U)*U;
% else if(T(j)>=tm/4-epsilon && T(j)<=tm/4+epsilon)
%     a = .1*a;
% b = .9*b;
%     RHSU = (SPMM-dt*gam*SPMM+dt*gam*SPC)*U+dt*gam*a*A;
%     RHSV = (SPMM-gam*dt*SPD)*V+dt*gam*b*A;
%     U =TMatrixU\RHSU;
%     V = TMatrixV\RHSV;
%     MatrixU(j,:)=U;
%     MatrixV(j,:)= V;
%     V = 150*V;
%     U = 600*U;    
% else if(T(j)>=tm/2-epsilon && T(j)<=tm/2+epsilon)
%         a = .5*T(j)*a;
%         b = .8*T(j)*b;
%     RHSU = (SPMM-dt*gam*SPMM+dt*gam*SPC)*U+dt*gam*a*A;
%     RHSV = (SPMM-gam*dt*SPD)*V+dt*gam*b*A;
%     U =TMatrixU\RHSU;
%     V = TMatrixV\RHSV;
%     MatrixU(j,:)=U;
%     MatrixV(j,:)= V;
%     V = 70*V;
%     U = 45*U;
%     else if(T(j)>=3*tm/4-epsilon && T(j)<=3*tm/4+epsilon)
%         a = .1*a;
%         b = .9*b;
%     RHSU = (SPMM-dt*gam*SPMM+dt*gam*SPC)*U+dt*gam*a*A;
%     RHSV = (SPMM-gam*dt*SPD)*V+dt*gam*b*A;
%     U =TMatrixU\RHSU;
%     V = TMatrixV\RHSV;
%     MatrixU(j,:)=U;
%     MatrixV(j,:)= V;
%     U = 23*U;
%      V = 23*V;
%         end
%     end
%     end
    %figure(1)
    
%subplot(1,2,1)
trisurf(LNODES,x,y,U(:,:))
 colorbar
shading interp
xlabel('x','fontsize',16) 
 view(2)
ylabel('y','fontsize',16)
zlabel('u','fontsize',16)
% title(['Evolution of u at t= ',num2str(T(j))],'fontsize',8)
title('Evolution of u at t = 9.0','fontsize',18)
axis equal tight
% subplot(1,2,2)
% trisurf(LNODES,x,y,1.5/max(U)*U(:,:))
% colorbar
% shading interp
% xlabel('x','fontsize',16) 
% grid off
% view(2)
% ylabel('y','fontsize',16)
% zlabel('u & v','fontsize',16)
% title(['Pattern formed by u at t = 30'])%at t= ',num2str(T(j))],'fontsize',8)
% axis equal tight
% MV(j)=getframe(gcf);
pause(1e-10) 


% end
%     RHSU = (SPMM-dt*gam*SPMM+dt*gam*SPC)*U+dt*gam*a*A;
%     RHSV = (SPMM-gam*dt*SPD)*V+dt*gam*b*A;
%     U = TMatrixU\RHSU;
%     V = TMatrixV\RHSV;
%     MatrixU(j,:)=U;
%     MatrixV(j,:)=V;
end
% 
% set(findobj('type','legend'),'fontsize',16)
% set(findobj('type','axes'),'fontsize',18)
% L2U = zeros(length(T),1);
% L2V = zeros(length(T),1);
% 
% for k = 1: length(T)-1
%    if (k>=1 && k<=60)
%     L2U(k) = sqrt(sum(abs(MatrixU(k+1,:)-MatrixU(k,:)).^2)/(T(j)-T(j-1)));
%     L2V(k) = sqrt(sum(abs(MatrixV(k+1,:)-MatrixV(k,:)).^2)/(T(j)-T(j-1)));
%    end
%     L2U(k) = sqrt(sum(abs(MatrixU(k+1,:)-MatrixU(k,:)).^2)/(T(j)-T(j-1)));
%     L2V(k) = sqrt(sum(abs(MatrixV(k+1,:)-MatrixV(k,:)).^2)/(T(j)-T(j-1)));
% end
% 
% L2U(length(T))=L2U(length(T)-1);
% L2V(length(T))=L2V(length(T)-1);
% figure(1)
% plot(T,.8/max(L2U)*L2U,'LineWidth',2,'color','r')
% 
% hold on
% plot(T,.75/max(L2V)*L2V,'LineWidth',2,'color','b')
% ylim([0 1.4]);
% set(findobj('type','legend'),'fontsize',20)
% set(findobj('type','axes'),'fontsize',20)
% legend('(||U^{m+1}-U^m||/\tau)','(||V^{m+1}-V^m||/\tau)','Location','NorthEast')
% xlabel('Time','fontsize',20)
% % ylabel('log(||\cdot||_L_2)','fontsize',20)
% ylabel('||\cdot||_L_2','fontsize',20)
% title('Convergence of solutions','fontsize',16)

% movie2avi(MV,'Shell_pattern.avi');
  set(findobj('type','legend'),'fontsize',18)
set(findobj('type','axes'),'fontsize',18)

