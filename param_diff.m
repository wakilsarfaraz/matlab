% This script simulates the parameter space for reaction system with
% diffusion with Schnakneberg Model.
% Author Wakil Sarfaraz 25/04/2016
%clear all;
tic
xmax =1;
N = 40; %(Number of points on the x, y interval on which the equation is solved.
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
%L = 100;
L = 15;
% T = (-(x+y).^3-x+y)./(x+y)-((d+1)*(n^2+m^2)*pi^2)/xmax^2;
% D = ((y-x)./(y+x)-((n^2+m^2)*pi^2)/xmax^2).*(-(y+x).^2-(d+1)*((n^2+m^2)*pi^2)/xmax^2)+2*y.*(y+x);
% Dcr = T.^2-4*D;
% img = sqrt(Dcr);




% 
% M = 100;
% dmax = 500;
% Diff = linspace(1,dmax,M+1);
%     for k = 1: length(Diff)
%  
gam = 1.;
% d = Diff(k);
d = 5;
T = gam*(-(x+y).^3-x+y)./(x+y)-((d+1)*(n^2+m^2)*pi^2)/L^2;
D = (gam*(y-x)./(y+x)-((n^2+m^2)*pi^2)/L^2).*(-gam*(y+x).^2-(d+1)*((n^2+m^2)*pi^2)/L^2)+2*gam^2*y.*(y+x);
Dcr = T.^2-4*D;
img = sqrt(Dcr);

Ru = zeros(NNODES, 2);
for i = 1: NNODES
    for j = 1 : 3
    %if (Dcr(LNODES(i,j))<=1e-2 && Dcr(LNODES(i,j))>=-1e-2 )
    %if(Dcr(LNODES(i,j))<0 && T(LNODES(i,j))<0)
    %if(Dcr(LNODES(i,j))>0 && T(LNODES(i,j))+sqrt(Dcr(LNODES(i,j)))>0)
    %if(Dcr(LNODES(i,j))<0 && T(LNODES(i,j))>0)
    if(T(LNODES(i,j))<=1e-3 && T(LNODES(i,j))>=-1e-3 && D(LNODES(i,j))>0)
        Ru(i,1) = x(LNODES(i,j));
        Ru(i,2) = y(LNODES(i,j));
    end
  
    end
    
end
plot(Ru(:,1),Ru(:,2),'.','Color','k')
%  hold on
%  text(0.4,1.8,'d=1.5','fontsize',15)
%  text(0.45,1.4,'d=2','fontsize',15)
%   text(0.2,1.05,'d=2.3','fontsize',15)
%   text(0.45,0.55,'d=2.8','fontsize',15)
%  text(0.45,0.22,'d=3.5','fontsize',15)

%  text(1.1,0.85,'d=2','fontsize',15)
%  text(0.95,0.85,'d=4','fontsize',15)
%   text(0.85,0.85,'d=6','fontsize',15)
%   text(0.75,0.85,'d=8','fontsize',15)
%  text(0.6,0.85,'d=10','fontsize',15)

%  text(0.6,0.1,'E','fontsize',20)
% text(0.55,0.45,'D','fontsize',20)
% text(0.35,0.85,'C','fontsize',20)
% text(0.1,1.4,'B','fontsize',20)
% text(0.15,1.75,'A','fontsize',20)

%  text(1,1.5,'A','fontsize',20)
% text(0.55,0.45,'E','fontsize',20)
% text(0.35,0.85,'D','fontsize',20)
% text(0.1,1.4,'C','fontsize',20)
% text(0.15,1.75,'B','fontsize',20)

%  text(0.0245,0.34,'A','fontsize',20)
% text(0.0225,0.33,'B','fontsize',20)
% text(0.0215,0.295,'C','fontsize',20)
% text(0.02,0.28,'D','fontsize',20)
% text(0.01,0.3,'E','fontsize',20)

% text(1.25,0.2,'E','fontsize',20)
% text(1.17,0.2,'D','fontsize',20)
% text(1.08,0.2,'C','fontsize',20)
% text(0.98,0.2,'B','fontsize',20)
% text(0.82,0.2,'A','fontsize',20)
% 
% text(1.25,0.2,'B','fontsize',20)
% text(1.17,0.2,'C','fontsize',20)
% text(1.08,0.2,'D','fontsize',20)
% text(0.98,0.2,'E','fontsize',20)
% text(1.34,0.2,'A','fontsize',20)

% text(0.025,0.35,'E','fontsize',20)
% text(0.07,0.35,'D','fontsize',20)
% text(0.105,0.35,'C','fontsize',20)
% text(0.13,0.35,'B','fontsize',20)
% text(0.147,0.35,'A','fontsize',20)

% text(0.001,0.67,'E','fontsize',20)
% text(0.025,0.59,'D','fontsize',20)
% text(0.06,0.49,'C','fontsize',20)
% text(0.085,0.4,'B','fontsize',20)
% text(0.105,0.4,'A','fontsize',20)

text(0.018,0.58,'c1','fontsize',20)
text(0.056,0.46,'c2','fontsize',20)
text(0.076,0.38,'c3','fontsize',20)
text(0.084,0.25,'c4','fontsize',20)
text(0.07,0.14,'c5','fontsize',20)

% % legend('d=1.5')
title('Transcritical bifurcation curves','fontsize',24)
%legend('Black:     d=2','Magenta:d=4','Green:    d=6','Red:       d=8','Blue:      d=10')
% title(['Unstable region with d=',num2str(Diff(k))],'fontsize',6)
xlabel('\alpha','fontsize',25)
ylabel('\beta','fontsize',25)
set(findobj('type','legend'),'fontsize',20)
set(findobj('type','axes'),'fontsize',20)
% 
%   xlim([0 1.4])
%    ylim([0 2.5])
 hold on

% Ru1 = zeros(NNODES, 2);
% Ru1 = zeros(NNODES, 2);
% for i = 1: NNODES
%     for j = 1 : 3
%     if (Dcr(LNODES(i,j))<=1e-3 && Dcr(LNODES(i,j))>=-1e-3 )
%         Ru1(i,1) = x(LNODES(i,j));
%         Ru1(i,2) = y(LNODES(i,j));
%     end
%   
%     end
%     
% end
% % 
% 
% % 
% figure(1)
% %subplot(2,2,1)
% plot(Ru1(:,1),Ru1(:,2),'.','Color','y')
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