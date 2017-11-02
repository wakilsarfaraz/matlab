% This script simulates the parameter space for reaction system without
% diffusion with Schnakneberg Model.
% Author Wakil Sarfaraz 25/04/2016
clear all;
tic
xmax =1;
N = 300; %(Number of points on the x, y interval on which the equation is solved.
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

rho = 35;
n = 10.3;
k = 0;
d =2;
gam = 1;
Trace = (-(x+y).^3-x+y)./(x+y)-((d+1)*(4*(2*k+1)*(n+2*k+1)*(n+4*k))/(rho^2*(n+4*k+2)));
Deter = (gam*(y-x)./(y+x)-(4*(2*k+1)*(n+2*k+1)*(n+4*k))/(rho^2*(n+4*k+2))).*...
    (-gam*(y+x).^2-(d+1)*(4*(2*k+1)*(n+2*k+1)*(n+4*k))/(rho^2*(n+4*k+2)))+2*gam^2*y.*(y+x);
Discrim = Trace.^2-4*Deter;
pureimag = sqrt(Discrim);

First = zeros(NNODES, 4);
for i = 1: NNODES
    for j = 1 : 3
        %if (Discrim(LNODES(i,j))< 0 && Trace(LNODES(i,j))>0)
    %if (Discrim(LNODES(i,j))< 1e-2 && Discrim(LNODES(i,j))>-1e-2)% && Trace(LNODES(i,j))<=0)
    if (Discrim(LNODES(i,j))>=0 && Trace(LNODES(i,j))+sqrt(Discrim(LNODES(i,j)))>=0)% && Trace(LNODES(i,j))-sqrt(Discrim(LNODES(i,j)))<0)
        First(i,1) = x(LNODES(i,j));
        First(i,2) = y(LNODES(i,j));
%     elseif (Discrim(LNODES(i,j))<0 && Trace(LNODES(i,j))<=1e-3 && Trace(LNODES(i,j))>=-1e-3)
%         First(i,3) = x(LNODES(i,j));
%         First(i,4) = y(LNODES(i,j));
    end
  
    end
    
end
plot(First(:,1),First(:,2),'.','Color','y')
%xlim([0 xmax/2])
% ylim([0.01 xmax])
title('Real unstable (Turing) region','fontsize',22)
xlabel('\alpha','fontsize',22)
ylabel('\beta','fontsize',22)
 hold on
%  plot(First(:,3),First(:,4),'.','Color','r')
%  legend('curve c_i','curve t_i')


 %%Real and Complex Partition Labels A
%   text(0.55,1.5,'d=2','fontsize',16)
%   text(0.32,1.5,'d=2.5','fontsize',16)
%   text(0.15,1.45,'d=3','fontsize',16)
%   text(0.25,1.15,'d=3.5','fontsize',16)
%    text(0.55,0.6,'d=4','fontsize',16)

 %%Real and Complex Partition Labels B
%   text(0.65,1.5,'d=2','fontsize',16)
%   text(0.74,1.2,'d=7','fontsize',16)
%   text(0.85,0.9,'d=12','fontsize',16)
%   text(0.95,0.6,'d=17','fontsize',16)
%    text(0.95,0.3,'d=22','fontsize',16)
 
%Real Region labels A:
%  text(0.46,0.66,'A','fontsize',16)
%   text(0.3,1.3,'B','fontsize',16)
%   text(0.25,1.5,'C','fontsize',16)
%   text(0.42,1.5,'D','fontsize',16)
%    text(1,2,'E','fontsize',16)

%%Real Region labels B:
%  text(1.4,0.25,'A','fontsize',16)
%   text(1.26,0.25,'B','fontsize',16)
%   text(1.19,0.25,'C','fontsize',16)
%   text(1.11,0.25,'D','fontsize',16)
%    text(1.02,0.25,'E','fontsize',16)

%Real Stable Region labels B:
%  text(1.38,0.2,'A','fontsize',16)
%   text(1.26,0.25,'B','fontsize',16)
%   text(1.19,0.25,'C','fontsize',16)
%   text(1.11,0.25,'D','fontsize',16)
%    text(1.02,0.25,'E','fontsize',16)

%Real Unstable (Turing) Region labels B:
%  text(0.14,0.32,'A','fontsize',16)
%   text(0.115,0.32,'B','fontsize',16)
%   text(0.085,0.32,'C','fontsize',16)
%   text(0.055,0.32,'D','fontsize',16)
%    text(0.025,0.2,'E','fontsize',16)


%%Complex Region labels A:
%  text(0.5,0.2,'A','fontsize',16)
%   text(0.3,0.9,'B','fontsize',16)
%   text(0.3,1.5,'D','fontsize',16)
%   text(0.42,1.55,'E','fontsize',16)
%    text(0.15,1.25,'C','fontsize',16)

%%Complex Unstable Region labels B:
%  text(0.007,0.67,'A','fontsize',16)
%   text(0.07,0.48,'B','fontsize',16)
%   text(0.097,0.43,'C','fontsize',16)
%   text(0.125,0.38,'D','fontsize',16)
%    text(0.147,0.38,'E','fontsize',16)

%%Complex Region labels B:
%  text(1.26,0.25,'A','fontsize',16)
%   text(1.19,0.25,'B','fontsize',16)
%   text(1.11,0.25,'C','fontsize',16)
%   text(1.02,0.25,'D','fontsize',16)
%    text(0.88,0.25,'E','fontsize',16)

%  text(0.05,0.5,'A','fontsize',16)
%   text(0.075,0.48,'B','fontsize',16)
%   text(0.085,0.48,'C','fontsize',16)
%   text(0.095,0.48,'D','fontsize',16)
%   text(0.105,0.48,'E','fontsize',16)

% 
% text(0.076,0.17,'c_1','fontsize',12)
%     text(0.085,0.23,'c_2','fontsize',12)
%     text(0.086,0.3,'c_3','fontsize',12)
%   text(0.082,0.36,'c_4','fontsize',12)
%     text(0.064,0.46,'c_5','fontsize',12)

% text(0.076,0.17,'c_1','fontsize',12)
%     text(0.085,0.23,'c_2','fontsize',12)
%     text(0.086,0.3,'c_3','fontsize',12)
%   text(0.082,0.36,'c_4','fontsize',12)
%     text(0.064,0.46,'c_5','fontsize',12)

%  text(0.029,0.55,'1','fontsize',16)
%      text(0.062,0.45,'2','fontsize',16)
%      text(0.082,0.35,'3','fontsize',16)
%    text(0.086,0.28,'4','fontsize',16)
%      text(0.086,0.23,'5','fontsize',16)
%       text(0.078,0.18,'6','fontsize',16)
%  text(0.07,0.15,'7','fontsize',16)
%   text(0.05,0.1,'8','fontsize',16) 
%   set(findobj('type','legend'),'fontsize',8)
 set(findobj('type','axes'),'fontsize',20)
%  
% hold on
%  