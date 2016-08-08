% This script simulates the parameter space for reaction system without
% diffusion with Schnakneberg Model.
% Author Wakil Sarfaraz 25/04/2016
clear all;
tic
xmax =3;
N = 5000; %(Number of points on the x, y interval on which the equation is solved.
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


m = 1;
n = 1;
L=21;
% T = (-(x+y).^3-x+y)./(x+y)-((d+1)*(n^2+m^2)*pi^2)/xmax^2;
% D = ((y-x)./(y+x)-((n^2+m^2)*pi^2)/xmax^2).*(-(y+x).^2-(d+1)*((n^2+m^2)*pi^2)/xmax^2)+2*y.*(y+x);
% Dcr = T.^2-4*D;
% img = sqrt(Dcr);




% 
% M = 100;
% dmax = 500;
% Diff = linspace(1,dmax,M+1);
%     for k = 1: length(Diff)
%         d = Diff(k);
d =27.5;
T = (-(x+y).^3-x+y)./(x+y)-((d+1)*(n^2+m^2)*pi^2)/L^2;
D = ((y-x)./(y+x)-((n^2+m^2)*pi^2)/L^2).*(-(y+x).^2-(d+1)*((n^2+m^2)*pi^2)/L^2)+2*y.*(y+x);
Dcr = T.^2-4*D;
img = sqrt(Dcr);

Ru = zeros(NNODES, 2);
for i = 1: NNODES
    for j = 1 : 3
%     if (T(LNODES(i,j))<=1e-3 && T(LNODES(i,j))>=-1e-3)
if (T(LNODES(i,j))<=1e-3 && T(LNODES(i,j))>=-1e-3 && Dcr(LNODES(i,j))<0)
        Ru(i,1) = x(LNODES(i,j));
        Ru(i,2) = y(LNODES(i,j));
    end
  
    end
    
end
% 
% 
% % 
% figure(1)
%subplot(2,2,1)
 plot(Ru(:,1),Ru(:,2),'.','Color','k')
 hold on
% text(0.0001,0.2,'Unstable in absence of diffusion, d=0','fontsize',10)
%  text(0.03,0.55,'c1','fontsize',16)
%  text(0.085,0.35,'c2','fontsize',16)
%  text(0.085,0.2,'c3','fontsize',16)
%  text(0.07,0.13,'c4','fontsize',16)
%  text(0.03,0.05,'c5','fontsize',16)
% text(0.06,0.55,'\leftarrow F','fontsize',13)
% text(0.001,0.675,'G','fontsize',13)
% % legend('d=1.5')
title('Boundaries for complex \lambda','fontsize',18)
% title(['Unstable region with d=',num2str(Diff(k))],'fontsize',6)
xlabel('\alpha','fontsize',18)
ylabel('\beta','fontsize',18)
set(findobj('type','legend'),'fontsize',8)
set(findobj('type','axes'),'fontsize',18)
% 
% xlim([0 0.15])
%  ylim([0 0.7])
% hold on

Ru1 = zeros(NNODES, 2);
Ru1 = zeros(NNODES, 2);
for i = 1: NNODES
    for j = 1 : 3
    if (Dcr(LNODES(i,j))<=1e-3 && Dcr(LNODES(i,j))>=-1e-3 )
        Ru1(i,1) = x(LNODES(i,j));
        Ru1(i,2) = y(LNODES(i,j));
    end
  
    end
    
end
% % 
% 
% % 
% figure(1)
% %subplot(2,2,1)
plot(Ru1(:,1),Ru1(:,2),'.','Color','y')
% % text(1,.1,'c1 \rightarrow','fontsize',16)
% %legend('Stable Spiral')
% title('Boundary Curves for complex \lambda','fontsize',16)
% xlabel('\alpha','fontsize',18)
% ylabel('\beta','fontsize',18)
% %  xlim([0 xmax/2])
% %  ylim([0 xmax])
% hold on
% pause(1e-10)
%     end
%   trisurf(LNODES,x,y,T(:,:))
%   title('The real part of complex \sigma','fontsize',18)  
%   shading interp
%   xlabel('\alpha','fontsize',16)
%   ylabel('\beta','fontsize',16)
%   zlabel('T(\alpha,\beta)','fontsize',16)
% %   text(0.014,0.1,'d=15','fontsize',12)
% 
%   set(findobj('type','legend'),'fontsize',8)
% set(findobj('type','axes'),'fontsize',12)
%  
% hold on
%  