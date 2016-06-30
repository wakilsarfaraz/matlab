% % This script simulates the parameter space for reaction system without
% % diffusion with Schnakneberg Model.
% % Author Wakil Sarfaraz 25/04/2016
% clear all;
% tic
% xmax =3;
% N = 3000; %(Number of points on the x, y interval on which the equation is solved.
% X = linspace(0,xmax,N+1); %This divides the interval into N equispaced sub intervals.
% %X = 0: 1/N :1; This can also be used to create the same X.
% [x, y] = meshgrid(X,X); % This creates an (N+1) by (N+1) grid of values ...
%                         % of X and assignes them to x and y. (So x and y
%                         % are now matrices of size (N+1)x(N+1)
% x =x(:);  % This turns the matrix x into a column vector.
% y =y(:);  % This turns the matrix y into a column vector.
% NNODES = (N+1)^2; % This introduces the number of nodes in the result of triangulation.
% NTRI = 2*N^2;  % This is the number of triangles in the mesh.
% LNODES = zeros(NTRI, 3); % This creates a zero matrix of the size 'number of triangles by 3'
%                          % (since there are 3 vertices for each triangle.
% 
% for i = 1:N
%     for j = 1: N
%         LNODES(i+2*(j-1)*N,1) = i+(j-1)*(N+1);%This numbers each of the first nodes of all the lower triangles in each square.
%         LNODES(i+2*(j-1)*N,2) = i+j*(N+1); % This numbers each of the second nodes of all the lower triangles in each square.
%         LNODES(i+2*(j-1)*N,3) = (i+1)+(j-1)*(N+1); % This numbers each of the third nodes of all the lower triangles in each square.
%         
%         LNODES(i+N+2*(j-1)*N,1) = i+1+j*(N+1); %This numbers each of the first nodes of all the upper triangles in each square.
%         LNODES(i+N+2*(j-1)*N,2) = (i+1)+(j-1)*(N+1);%This numbers each of the second nodes of all the upper triangles in each square.
%         LNODES(i+N+2*(j-1)*N,3) = i+j*(N+1);%This numbers each of the third nodes of all the upper triangles in each square.
%                                             % Here for those nodes that
%                                             % some triangles share (copy
%                                             % and paste can work).
%     end
% end
% 
% T = (-(x+y).^3-x+y)./(x+y);
% D = (x+y).^2;
% C = y-x+(y+x).^3-2*(y+x).^2;
% Dcr = T.^2-4*D;
% 
% Ru = zeros(NNODES, 2);
% for i = 1: NNODES
%     for j = 1 : 3
%          if (T(LNODES(i,j))<=-1e-4 )%&& T(LNODES(i,j))<=-1e-2)
% %     if (Dcr(LNODES(i,j))>0 && sqrt(Dcr(LNODES(i,j)))+T(LNODES(i,j)) < 0)
%         Ru(i,1) = x(LNODES(i,j));
%         Ru(i,2) = y(LNODES(i,j));
%     end
%   
%     end
%     
% end
% 
% 
% % 
% figure(1)
% %subplot(2,2,1)
% plot(Ru(:,1),Ru(:,2),'.','Color','b')
% %legend('Stable Spiral')
% title('Stability partition','fontsize',16)
% xlabel('\alpha','fontsize',18)
% ylabel('\beta','fontsize',18)
%  xlim([0 0.3])
%  ylim([0 1.2])
% hold on
% 
% Ru1 = zeros(NNODES, 2);
% for i = 1: NNODES
%     for j = 1 : 3
%     if (T(LNODES(i,j))>=1e-4)%&& T(LNODES(i,j))<=-1e-2)
% %     if (Dcr(LNODES(i,j))>0 && sqrt(Dcr(LNODES(i,j)))+T(LNODES(i,j)) < 0)
%         Ru1(i,1) = x(LNODES(i,j));
%         Ru1(i,2) = y(LNODES(i,j));
%     end
%   
%     end
%     
% end
% plot(Ru1(:,1),Ru1(:,2),'.','Color','r')
% 
% 
% Ru2 = zeros(NNODES, 2);
% for i = 1: NNODES
%     for j = 1 : 3
%      if (T(LNODES(i,j))<1e-2 && T(LNODES(i,j))>-1e-2)
% %     if (Dcr(LNODES(i,j))>0 && sqrt(Dcr(LNODES(i,j)))+T(LNODES(i,j)) < 0)
%         Ru2(i,1) = x(LNODES(i,j));
%         Ru2(i,2) = y(LNODES(i,j));
%     end
%   
%     end
%     
% end
% plot(Ru2(:,1),Ru2(:,2),'.','Color','y')
% legend('Blue = Stable region','Red = Unstable region','Yellow = Centre (Partition)')
% 
x = [0.02 0.051 0.096 0.1 0.6 1 1.4];
y = [0.13 0.105 0.102 0.5 1 1 2];

plot(x,y,'*')
grid on
% xlim([-0.3 2])
% ylim([0 3])
% axis tight
X = x+y;
Y = y./(x+y).^2;

plot(X,Y,'*')
grid on