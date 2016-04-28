% This script simulates the parameter space for reaction diffusion system
% with Schnakneberg Model.
% Author Wakil Sarfaraz 27/04/2016
clear all;
tic
xmax =1;
N = 100; %(Number of points on the x, y interval on which the equation is solved.
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


M = 100;
dmax = 100;
D = linspace(0,dmax,M-1);

for n = 1: 1
    for k = 1: length(D)
        d = D(k)
A = (x+y).^2-2*(x.^2-y.^2);
limitu =[min(A) max(A)]
Ru = zeros(NNODES, 2);
for i = 1: NNODES
    for j = 1 : 3
    if (A(LNODES(i,j)) > 0)
        Ru(i,1) = x(LNODES(i,j));
        Ru(i,2) = y(LNODES(i,j));
    end
  
    end
    

end
figure(1)
subplot(2,2,1)
plot(Ru(:,1),Ru(:,2),'.','Color','g')
title(['A is positive with n=',num2str(n),' and d=',num2str(D(k))],'fontsize',6)
xlabel('alpha')
ylabel('beta')

B = 2*n^2*pi^2*(d*(y-x)-(x+y).^3)./(y+x);
 limitg =[min(B) max(B)]
Rg1 = zeros(NNODES, 2);
for i = 1: NNODES
    for j = 1 : 3
    if (B(LNODES(i,j)) > 0)
        Rg1(i,1) = x(LNODES(i,j));
        Rg1(i,2) = y(LNODES(i,j));
    end
  
    end
    
end
%figure(2)
subplot(2,2,2)
plot(Rg1(:,1),Rg1(:,2),'.','Color','r')
title(['Necessary for instability with n=',num2str(n),' and d=',num2str(D(k))],'fontsize',6 )
xlabel('alpha')
ylabel('beta')


Rg2 = zeros(NNODES,2);
for i = 1 : NNODES
    for j = 1 : 3
        if (abs(B(LNODES(i,j))) > abs(A(LNODES(i,j)))+4*n^4*pi^4)
            Rg2(i,1) = x(LNODES(i,j));
            Rg2(i,2) = y(LNODES(i,j));
        end
    end
end
%figure(3)
subplot(2,2,3)
plot(Rg2(:,1),Rg2(:,2),'.','Color','b')
title(['Sufficient for instability with n=',num2str(n),' and d=',num2str(D(k))],'fontsize',6 )
xlabel('alpha')
ylabel('beta')

subplot(2,2,4)
trisurf(LNODES,x,y,A)
shading interp
title('B plane','fontsize',6)

pause(1e-1)
    end

    end
