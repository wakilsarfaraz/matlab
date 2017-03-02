clear all; close all; clc;
%% This program solves Poisson's equation on a disk.
% equation is of the form laplace u = -pi^2 (x^2+y^2)cos(pi/2*(x^2+y^2))
% Wakil Sarfaraz 02/02/2017
clear all;
tic
xmax = 1;
N = 10;

rho = 1;
nu = 0.3;
k=1;

 fd=inline('sqrt(sum(p.^2,2))-0.5','p');
 [p,t]=distmesh2d(fd,@huniform,xmax/N,[-1,-1;1,1],[]);

x = p(:,1);
y = p(:,2);

NNODES = length(x);
NTRI = size(t,1);
LNODES = t;


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
       
       F(1) = (1-ksi-eta)*(pi^2*(xx.^2+yy.^2)*cos(pi/2*(xx.^2+yy.^2))+2*pi*sin(pi/2*(xx.^2+yy.^2)))*det(J)*1/2; %This computes the integral for the first, second and third                                                              %vertix of traingles
       F(2) = (ksi)*(pi^2*(xx.^2+yy.^2)*cos(pi/2*(xx.^2+yy.^2))+2*pi*sin(pi/2*(xx.^2+yy.^2)))*det(J)*1/2;       % In the transformed coordinates.
       F(3) = (eta)*(pi^2*(xx.^2+yy.^2)*cos(pi/2*(xx.^2+yy.^2))+2*pi*sin(pi/2*(xx.^2+yy.^2)))*det(J)*1/2;


       for i = 1 : 3
         
           LV(LNODES(n,i)) = LV(LNODES(n,i))+ F(i); % This assigns all the newly computed values to the Load vector.

      
       end 
     
end
for i = 1: NNODES
  if (abs(fd(p(i,:)))<=1e-5 )

        LV(i) = 1;
        SP(i,:) = 0;
        SP(i,i) =1 ;

    end
end
U = SP\LV; %Solves the linear system.

toc

 u = cos(pi/2*(x.^2+y.^2));
 u = u-min(u);

 


figure(1)

subplot(3,1,1)
%figure(1)
trisurf(LNODES,x,y,u)
%view(2)
xlim([-1 1]);
ylim([-1 1]);
%zlim([0.925,1]);
xlabel('x') 
ylabel('y')
zlabel('u(x,y)')
title('Exact Solution')
shading interp

 subplot(3,1,2)
% figure(2)
trisurf(LNODES,x,y,U)
%view(2)
xlim([-1 1]);
ylim([-1 1]);
%zlim([0 1]);
xlabel('x') 
ylabel('y')
zlabel('U(x,y)')
title('Numerical Solution')
shading interp
D = abs(U-u);
C = u-U;
%figure(3)
subplot(3,1,3)
trisurf(LNODES,x,y,C)
title('Error')
% 
 Error = sum(D.^2)
%  
%  
%  figure(2)
% trisurf(LNODES,x,y,D)
% shading interp
