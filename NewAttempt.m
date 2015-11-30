% This program solves Poisson's equation in Cartesian coordinates.
% Wakil Sarfaraz with Help of Dr Kathryn Gillow. 22 May 2014 
clear all;
tic
xmax =1;
N = 50; %(Number of points on the x, y interval on which the equation is solved.
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

SP = sparse(NNODES, NNODES);  %This creates a sparse matrix of the size 'Number of nodes by Number of nodes'(entries usually computed by integration).
LV = zeros(NNODES,1);         % This creates the load vector of zero entries for now which is l(v) part of the weak formulation.

for n = 1: NTRI
    r1 = [x(LNODES(n,1)) y(LNODES(n,1))];% Position vector (x,y)' for the first nodes of all triangles.
    r2 = [x(LNODES(n,2)) y(LNODES(n,2))];% Position vector (x,y)' for the second nodes of all triangles.
    r3 = [x(LNODES(n,3)) y(LNODES(n,3))];% Position vector (x,y)' for the third nodes of all triangles.
    J = [r2(1)-r1(1) r2(2)-r1(2); r3(1)-r1(1) r3(2)-r1(2)]; %This is the jacobian matrix of the mapping.
    Astiff = (1/(2*det(J)))* [(r2-r3)*(r2-r3)' (r2-r3)*(r3-r1)' (r2-r3)*(r1-r2)';... %This is the Global Stiffness Matrix (Lecture Notes)
           (r2-r3)*(r3-r1)' (r3-r1)*(r3-r1)' (r3-r1)*(r1-r2)';...
           (r2-r3)*(r1-r2)' (r3-r1)*(r1-r2)' (r1-r2)*(r1-r2)'];
       for i = 1 : 3 
           for j = 1:3  % Since SP was created a zero matrix originally, Now for all three vertices of all triangles we put the 
               SP(LNODES(n,i),LNODES(n,j)) =SP(LNODES(n,i),LNODES(n,j))+ Astiff(i,j); %The values of Astiff in the SP.
           end
       end
       % Up to here we constructed the Matrix of the coefficients for the
       % Linear system. Now we start by computing the load vector.
       ksi =1/3; %This is the new variable of mapping. Value (1/3) is chosen due to quadrature rule.
       eta = 1/3; % This is the new variable of mapping. Value (1/3) is chosen due to quadrature rule.
       xx = (1-ksi-eta)*r1(1) + ksi* r2(1) + eta*r3(1); %This is x(ksi,eta) which is used for transformation.
       yy = (1-ksi-eta)*r1(2) + ksi* r2(2) + eta*r3(2); % This is y(ksi,eta) which is used for transformation.
       
%        F(1) = (1-ksi-eta)*5*pi^2*sin(pi*xx)*sin(2*pi*yy)*det(J)*1/2; %This computes the integral for the first, second and third
%                                                                       %vertix of traingles
%        F(2) = (ksi)*5*pi^2*sin(pi*xx)*sin(2*pi*yy)*det(J)*1/2;       % In the transformed coordinates.
%        F(3) = (eta)*5*pi^2*sin(pi*xx)*sin(2*pi*yy)*det(J)*1/2;
       
       F(1) = (1-ksi-eta)*det(J)*1/2;%| these will solve for laplace u =1 with zero dirichlet bcs.
       F(2) = ksi* det(J)*1/2;
       F(3) = eta* det(J)*1/2;
       for i = 1 : 3
           LV(LNODES(n,i)) = LV(LNODES(n,i))+ F(i); % This assigns all the newly computed values to the Load vector.
       end 
     
end
for i = 1: NNODES
    if (x(i)==0 || y(i)==0 || x(i)==xmax || y(i)==xmax) % This enforces boundary conditions.
        LV(i) = 0;
        SP(i,:) = 0;
        SP(i,i) =1 ;
    end
end

U = SP\LV; %Solves the linear system.
toc

u = sin(pi*x).*sin(2*pi*y);



subplot(1,2,1)
%figure(1)
trisurf(LNODES,x,y,u)
xlabel('x') 
ylabel('y')
zlabel('u(x,y)')
title('Exact Solution')
subplot(1,2,2)
%figure(2)
trisurf(LNODES,x,y,U)
xlabel('x') 
ylabel('y')
zlabel('u(x,y)')
title('Numerical Solution')
D = abs(U-u);


 Error = sum(D.^2)
    



  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
       
  
  
  
        