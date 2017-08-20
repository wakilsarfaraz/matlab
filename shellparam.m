% This script simulates the parameter space for reaction system without
% diffusion with Schnakneberg Model.
% Author Wakil Sarfaraz 25/04/2016
clear all;
tic
xmax =10;
N = 2000; %(Number of points on the x, y interval on which the equation is solved.
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

rho = 0.5;
a = 0.5;
l = .27;
k = 0;
d =3;
gam = 1;

Trace = (-(x+y).^3-x+y)./(x+y)-((d+1).*eta(a,rho,k,l).^2);
Deter = (gam*(y-x)./(y+x)-eta(a,rho,k,l).^2).*...
    (-gam*(y+x).^2-(d+1)*eta(a,rho,k,l).^2)+2*gam^2*y.*(y+x);
Discrim = Trace.^2-4*Deter;
pureimag = sqrt(Discrim);

First = zeros(NNODES, 4);

for i = 1: NNODES
    for j = 1 : 3
    if(Discrim(LNODES(i,j)) <= 1e-2 && Discrim(LNODES(i,j)) >= -1e-2)
        First(i,1) = x(LNODES(i,j));
        First(i,2) = y(LNODES(i,j));
    end
    end
end
plot(First(:,1),First(:,2),'.','Color','b')
 xlim([0 1.8])
 ylim([0 2.5])
title('Partition for region of real and complex \sigma_{1,2}','fontsize',14)
xlabel('\alpha','fontsize',20)
ylabel('\beta','fontsize',20)
hold on
%  hold on
%  plot(First(:,3),First(:,4),'.','Color','r')
 %legend('Partition in Blue: d=4','Partition in Red: d=2')
%  text(0.46,0.66,'A','fontsize',12)
%   text(0.3,1.3,'B','fontsize',12)
%   text(0.25,1.5,'C','fontsize',12)
%   text(0.42,1.5,'D','fontsize',12)
%    text(1,2,'F','fontsize',12)

%  text(0.5,0.2,'A','fontsize',12)
%   text(0.3,0.9,'B','fontsize',12)
%   text(0.3,1.5,'D','fontsize',12)
%   text(0.42,1.55,'F','fontsize',12)
%    text(0.15,1.25,'C','fontsize',12)

%  text(0.05,0.5,'A','fontsize',12)
%   text(0.075,0.48,'B','fontsize',12)
%   text(0.085,0.48,'C','fontsize',12)
%   text(0.095,0.48,'D','fontsize',12)
%   text(0.105,0.48,'F','fontsize',12)

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

% % 
  set(findobj('type','legend'),'fontsize',8)
 set(findobj('type','axes'),'fontsize',20)
%  
% hold on
%  