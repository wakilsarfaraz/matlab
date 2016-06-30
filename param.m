% This script simulates the parameter space for reaction system without
% diffusion with Schnakneberg Model.
% Author Wakil Sarfaraz 25/04/2016
clear all;
tic
xmax =1;
N = 1000; %(Number of points on the x, y interval on which the equation is solved.
X = linspace(0,xmax,N+1); %This divides the interval into N equispaced sub intervals.
%X = 0: 1/N :1; This can also be used to create the same X.
[x, y] = meshgrid(X,X); % This creates an (N+1) by (N+1) grid of values ...
                        % of X and assignes them to x and y. (So x and y
                        % are now matrices of size (N+1)x(N+1)
x =x(:);  % This turns the matrix x into a column vector.
y =y(:);  % This turns the matrix y into a column vector.
NNODES = (N+1)^2; % This introduces the number of nodes in the result of triangulation.
NTRI = 2*N^2;  % This is the number of triangles in the mesh.
LNODES = zeros(NTRI, 3); % This creates a zero matrix of the size 'number of triangles by 3'
                         % (since there are 3 vertices for each triangle.

for i = 1:N
    for j = 1: N
        LNODES(i+2*(j-1)*N,1) = i+(j-1)*(N+1);%This numbers each of the first nodes of all the lower triangles in each square.
        LNODES(i+2*(j-1)*N,2) = i+j*(N+1); % This numbers each of the second nodes of all the lower triangles in each square.
        LNODES(i+2*(j-1)*N,3) = (i+1)+(j-1)*(N+1); % This numbers each of the third nodes of all the lower triangles in each square.
        
        LNODES(i+N+2*(j-1)*N,1) = i+1+j*(N+1); %This numbers each of the first nodes of all the upper triangles in each square.
        LNODES(i+N+2*(j-1)*N,2) = (i+1)+(j-1)*(N+1);%This numbers each of the second nodes of all the upper triangles in each square.
        LNODES(i+N+2*(j-1)*N,3) = i+j*(N+1);%This numbers each of the third nodes of all the upper triangles in each square.
                                            % Here for those nodes that
                                            % some triangles share (copy
                                            % and paste can work).
    end
end


T = (-(x+y).^3-x+y)./(x+y);
D = (x+y).^2;
C = y-x+(y+x).^3-2*(y+x).^2;
Dcr = T.^2-4*D;

Tr = zeros(NNODES, 2);
for i = 1: NNODES
    for j = 1 : 3

        Tr(i,1) = T(LNODES(i,j));
        Tr(i,2) = Dcr(LNODES(i,j));
%     end
%   
    end
%     
end


% % 
% figure(1)
% %subplot(2,2,1)
plot(Tr(:,1),Tr(:,2),'.','Color','b')

% % text(.6,1,'A','fontsize',16)
% %legend('Stable Spiral')
% title('Stability plot for real \lambda','fontsize',16)
% xlabel('\alpha','fontsize',18)
% ylabel('\beta','fontsize',18)
%  xlim([0 1.36])
%  ylim([0 2.5])
% hold on

  

  
%   
% Ru1 = zeros(NNODES, 2);
% for i = 1: NNODES
%     for j = 1 : 3
%     % if (Dcr(LNODES(i,j))>=-1e-4)
% %      if (abs(T(LNODES(i,j)))<=1e-2)
%          if (Dcr(LNODES(i,j))>0 && sqrt(Dcr(LNODES(i,j)))+T(LNODES(i,j))>0)
%         Ru1(i,1) = x(LNODES(i,j));
%         Ru1(i,2) = y(LNODES(i,j));
%     end
%   
%     end
%     
% end
%    plot(Ru1(:,1),Ru1(:,2),'.','Color','y')
%    text(0.08,0.3,'\leftarrow c2','fontsize', 16)
%    text(.8,1,'c1 \rightarrow','fontsize',16)
%    text(1.2,1.5,'A','fontsize',16)
%    text(0,0.2,'B','fontsize',16)
%legend('Unstable spiral')
% Rg1 = zeros(NNODES, 2);
% for i = 1: NNODES
%     for j = 1 : 3
%         if( Dcr(LNODES(i,j)) <= 1e-2 && Dcr(LNODES(i,j))>= -1e-2 && T(LNODES(i,j))>0)
% % if (T(LNODES(i,j))>=-1e-2 && T(LNODES(i,j)) <= 1e-2)
%         Rg1(i,1) = x(LNODES(i,j));
%         Rg1(i,2) = y(LNODES(i,j));
%     end
%   
%     end
%     
% end
%figure(2)
% subplot(1,2,1)
% trisurf(LNODES,x,y,g)
% subplot(1,2,2)
% subplot(2,2,2)
%   plot(Rg1(:,1),Rg1(:,2),'.','Color','r')
  
%     Rg2 = zeros(NNODES, 2);
% for i = 1: NNODES
%     for j = 1 : 3
% %         if( Dcr(LNODES(i,j)) < 0)
% %     if (Dcr(LNODES(i,j))<= 1e-2 && Dcr(LNODES(i,j))>= -1e-2 )% && abs(T(LNODES(i,j)))< abs(sqrt(Dcr(LNODES(i,j)))))
%         if (Dcr(LNODES(i,j)) <= 1e-2 && Dcr(LNODES(i,j))>= -1e-2 && T(LNODES(i,j))<0)
%         Rg2(i,1) = x(LNODES(i,j));
%         Rg2(i,2) = y(LNODES(i,j));
%     end
%   
%     end
%     
% end
% %figure(2)
% % subplot(1,2,1)
% % trisurf(LNODES,x,y,g)
% % subplot(1,2,2)
% % subplot(2,2,2)
%   plot(Rg2(:,1),Rg2(:,2),'.','Color','b')
  

%title('Two distict real roots both negative','fontsize',6)
% xlabel('alpha')
% ylabel('beta')
% xlim([-0.1 0.6])
% ylim([-0.1 1.1])
%  
%  
%  
%  
%  
%  
%  Rg3 = zeros(NNODES,2);
%  for i = 1: NNODES
%     for j = 1 : 3
%     if (Dcr(LNODES(i,j))<0 && T(LNODES(i,j)) <0)
%         Rg3(i,1) = x(LNODES(i,j));
%         Rg3(i,2) = y(LNODES(i,j));
%     end
%   
%     end
%     
% end
% % %figure(3)
% % % subplot(1,2,1)
% % % trisurf(LNODES,x,y,g)
% % % subplot(1,2,2)
% % subplot(2,2,3)
%  plot(Rg3(:,1),Rg3(:,2),'.','Color','m')
%  text(0.6,1,'C','fontsize', 16)
%  legend('Blue = Unstable spiral','Megenta = Stable spiral','Red = Stable node','Green = Unstable node')
% legend('Stable node')
% title('Two distinct real roots with different signs','fontsize',6)
% xlabel('alpha')
% ylabel('beta')
% 
%  Rg4 = zeros(NNODES,2);
%  for i = 1: NNODES
%     for j = 1 : 3
%     if (Dcr(LNODES(i,j))<0 && T(LNODES(i,j)) >0)
%         Rg4(i,1) = x(LNODES(i,j));
%         Rg4(i,2) = y(LNODES(i,j));
%     end
%   
%     end
%     
% end
% % %figure(4)
% % % subplot(1,2,1)
% % % trisurf(LNODES,x,y,g)
% % % subplot(1,2,2)
% % subplot(2,2,4)
% plot(Rg4(:,1),Rg4(:,2),'.','Color','b')
% text(0.05,0.6,'D','fontsize', 16)
% title('Two repeated real e-values','fontsize',6)
% xlabel('alpha')
% ylabel('beta')

%  Rg5 = zeros(NNODES,2);
%  for i = 1: NNODES
%     for j = 1 : 3
%     if (T(LNODES(i,j))<0)
%         Rg5(i,1) = x(LNODES(i,j));
%         Rg5(i,2) = y(LNODES(i,j));
%     end
%   
%     end
%     
% end
% %figure(4)
% % subplot(1,2,1)
% % trisurf(LNODES,x,y,g)
% % subplot(1,2,2)
% subplot(2,2,4)
% plot(Rg5(:,1),Rg5(:,2),'.','Color','b')
% hold on
% title('Boundaries of complex region for \lambda','fontsize',14)
% text(0.8,1,'Stable spiral','fontsize', 12)
% hold on
% title('Two repeated real e-values','fontsize',6)
% xlabel('alpha')
% ylabel('beta')
%   Rg6 = zeros(NNODES,2);
%  for i = 1: NNODES
%     for j = 1 : 3
%     if (T(LNODES(i,j)) >0)
%         Rg6(i,1) = x(LNODES(i,j));
%         Rg6(i,2) = y(LNODES(i,j));
%     end
%   
%     end
%     
% end
% %figure(4)
% % subplot(1,2,1)
% % trisurf(LNODES,x,y,g)
% % subplot(1,2,2)
% subplot(2,2,4)
% plot(Rg6(:,1),Rg6(:,2),'.','Color','r')
% title('Stability partition','fontsize',16)
% xlabel('\alpha','fontsize',18)
% ylabel('\beta','fontsize',18)
%  xlim([0 0.25])
%  ylim([0 1.1])
% 
%    text(1.,1.5,' A','fontsize', 14)
%    text(0.01,0.2,'B','fontsize',14)
%    text(0.85,1,'c1 \rightarrow ','fontsize',12)
%    text(0.06,0.2,'\leftarrow c2','fontsize',12)
% % text(0.025,0.35,'\leftarrow c3 ','fontsize', 14)
% % title('Two repeated real e-values','fontsize',6)
% % xlabel('alpha')
% % ylabel('beta')
%   Rg7 = zeros(NNODES,2);
%  for i = 1: NNODES
%     for j = 1 : 3
%     if (Dcr(LNODES(i,j)) <= 1e-2 && Dcr(LNODES(i,j))>= -1e-2 && T(LNODES(i,j))>0)
%         Rg7(i,1) = x(LNODES(i,j));
%         Rg7(i,2) = y(LNODES(i,j));
%     end
%   
%     end
%     
% end
%figure(4)
% subplot(1,2,1)
% trisurf(LNODES,x,y,g)
% subplot(1,2,2)
% % subplot(2,2,4)
%  plot(Rg7(:,1),Rg7(:,2),'.','Color','r')
% 
% % text(.2,0.5,'\leftarrow c2','fontsize', 14)
% legend('A: \lambda_1 & \lambda_2 both real negative','B: \lambda_1 or \lambda_2 or both real positive',...
%     'c1: Stable star','c2: Unstable star')
% set(findobj('type','legend'),'fontsize',10)
% set(findobj('type','axes'),'fontsize',12)
% hold off
% title('Two repeated real e-values','fontsize',6)
% xlabel('alpha')
% ylabel('beta')
% trisurf(LNODES,x,y,T(:,:))
% shading interp
% max(T)