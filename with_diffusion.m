%% This script solves the time dependent problem using uniform triangulation.
% Wakil Sarfaraz   30/11/15
clear all;
tic
xmax = 1;
tm = 1;
dt = 0.01;
M = tm/dt;

N = 30; %(Number of points on the x, y interval on which the equation is solved.
X = linspace(0,xmax,N+1);
T = linspace(0, tm,M+1); %This divides the interval into N equispaced sub intervals.
%X = 0: 1/N :1; This can also be used to create the same X.
[x, y] = meshgrid(X,X); % This creates an (N+1) by (N+1) grid of values ...
                        % of X and assignes them to x and y. (So x and y
                        % are now matrices of size (N+1)x(N+1)
x =x(:);  % This turns the matrix x into a column vector.
y =y(:);  % This turns the matrix y into a column vector.
NNODES = (N+1)^2;
Q = ones(NNODES,1);% This introduces the number of nodes in the result of triangulation.
NTRI = 2*N^2;  % This is the number of triangles in the mesh.
LNODES = zeros(NTRI,3); % This creates a zero matrix of the size 'number of triangles by 3'
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

SP = sparse(NNODES, NNODES);
SPM = sparse(NNODES, NNODES);%This creates a sparse matrix of the size 'Number of nodes by Number of nodes'(entries usually computed by integration).
LV = zeros(NNODES,1);         % This creates the load vector of zero entries for now which is l(v) part of the weak formulation.

for n = 1: NTRI
    r1 = [x(LNODES(n,1)) y(LNODES(n,1))];% Position vector (x,y)' for the first nodes of all triangles.
    r2 = [x(LNODES(n,2)) y(LNODES(n,2))];% Position vector (x,y)' for the second nodes of all triangles.
    r3 = [x(LNODES(n,3)) y(LNODES(n,3))];% Position vector (x,y)' for the third nodes of all triangles.
    J = [r2(1)-r1(1) r2(2)-r1(2); r3(1)-r1(1) r3(2)-r1(2)]; %This is the jacobian matrix of the mapping.
%     ksi = 1/3;
%     eta = 1/3;
   Astiff = (((r2(1)+r3(1))/(2*det(J)))*...
       [(r2(2)-r3(2))^2+(r3(1)-r2(1))^2 ...
       (r2(2)-r3(2))*(r3(2)-r1(2))+(r3(1)-r2(1))*(r1(1)-r3(1))...
       (r2(2)-r3(2))*(r1(2)-r2(2))+(r3(1)-r2(1))*(r2(1)-r1(1)); ...
       (r3(2)-r1(2))*(r2(2)-r3(2))+(r1(1)-r3(1))*(r3(1)-r2(1)) ...
       (r3(2)-r1(2))^2+(r1(1)-r3(1))^2  ...
       (r3(2)-r1(2))*(r1(2)-r2(2))+(r1(1)-r3(1))*(r2(1)-r1(1)); ...
       (r2(2)-r3(2))*(r1(2)-r2(2))+(r3(1)-r2(1))*(r2(1)-r1(1))  ...
       (r3(2)-r1(2))*(r1(2)-r2(2))+(r1(1)-r3(1))*(r2(1)-r1(1))  ... 
       (r1(2)-r2(2))^2+(r2(1)-r1(1))^2]+(1/(2*det(J)))* [(r2-r3)*(r2-r3)' (r2-r3)*(r3-r1)' (r2-r3)*(r1-r2)';...
       (r2-r3)*(r3-r1)' (r3-r1)*(r3-r1)' (r3-r1)*(r1-r2)';...
       (r2-r3)*(r1-r2)' (r3-r1)*(r1-r2)' (r1-r2)*(r1-r2)']);
% Astiff = (1/(4*det(J)))* [(r2-r3)*(r2-r3)' (r2-r3)*(r3-r1)' (r2-r3)*(r1-r2)';... %This is the Global Stiffness Matrix (Lecture Notes)
%            (r2-r3)*(r3-r1)' (r3-r1)*(r3-r1)' (r3-r1)*(r1-r2)';...
%            (r2-r3)*(r1-r2)' (r3-r1)*(r1-r2)' (r1-r2)*(r1-r2)'];
% %    Astiff = (1/(4*det(J)))*...
%        [(r2(2)-r3(2))^2+(r3(1)-r2(1))^2 ...
%        (r2(2)-r3(2))*(r3(2)-r1(2))+(r3(1)-r2(1))*(r1(1)-r3(1))...
%        (r2(2)-r3(2))*(r1(2)-r2(2))+(r3(1)-r2(1))*(r2(1)-r1(1)); ...
%        (r3(2)-r1(2))*(r2(2)-r3(2))+(r1(1)-r3(1))*(r3(1)-r2(1)) ...
%        (r3(2)-r1(2))^2+(r1(1)-r3(1))^2  ...
%        (r3(2)-r1(2))*(r1(2)-r2(2))+(r1(1)-r3(1))*(r2(1)-r1(1)); ...
%        (r2(2)-r3(2))*(r1(2)-r2(2))+(r3(1)-r2(1))*(r2(1)-r1(1))  ...
%        (r3(2)-r1(2))*(r1(2)-r2(2))+(r1(1)-r3(1))*(r2(1)-r1(1))  ... 
%        (r1(2)-r2(2))^2+(r2(1)-r1(1))^2];
% %        
       
%    Amass = (det(J)/120)*[6*r1(1)+2*r2(1)+2*r3(1) 2*r1(1)+2*r2(1)+r3(1) 2*r1(1)+r2(1)+2*r3(1);...
%        2*r1(1)+2*r2(1)+r3(1) 2*r1(1)+2*r2(1)+2*r3(1) r1(1)+2*r2(1)+2*r3(1);...
%        2*r1(1)+r2(1)+2*r3(1) r1(1)+2*r2(1)+2*r3(1)  2*r1(1)+2*r2(1)+2*r3(1)];
   
      Amass = (det(J)/120)*[6*(r1+2*r2+2*r3)*(r1+2*r2+2*r3)' 2*(r1+2*r2+r3)*(r1+2*r2+r3)' (2*r1+r2+2*r3)*(2*r1+r2+2*r3)';...
       (2*r1+2*r2+r3)*(2*r1+2*r2+r3)' (2*r1+2*r2+2*r3)*(2*r1+2*r2+2*r3)' (r1+2*r2+2*r3)*(r1+2*r2+2*r3)';...
       (2*r1+r2+2*r3)*(2*r1+r2+2*r3)' (r1+2*r2+2*r3)*(r1+2*r2+2*r3)'  (2*r1+2*r2+2*r3)*(2*r1+2*r2+2*r3)'];
   

  
     

   
       for i = 1 : 3 
           for j = 1:3  % Since SP was created a zero matrix originally, Now for all three vertices of all triangles we put the 
               SP(LNODES(n,i),LNODES(n,j)) =SP(LNODES(n,i),LNODES(n,j))+ Astiff(i,j);
               SPM(LNODES(n,i),LNODES(n,j)) =SPM(LNODES(n,i),LNODES(n,j))+ Amass(i,j);%The values of Astiff in the SP.
           end
       end
       
    
       
end


TMatrix =  (dt*SP+SPM);



for i = 1: NNODES
    if (x(i)==xmax && y(i)>=0 && y(i) <= xmax)  
        LV(i) = 0;
        TMatrix(i,:) = 0;
        SPM(i,:) = 0;
        TMatrix(i,i) =1 ;
    elseif (x(i)==0 && y(i)>=0 && y(i) <= xmax)  
        LV(i) = 0;
        TMatrix(i,:) = 0;
        SPM(i,:) = 0;
        TMatrix(i,i) =1 ;
     elseif (y(i) == 0 && x(i) >= 0 && x(i) <= xmax) 
        TMatrix(i,:) = 0;
        TMatrix(i,i) = 1;
        SPM(i,:) = 0;
        LV(i) = 0;
    elseif ( y(i) == xmax && x(i) >= 0 && x(i) <= xmax) 
        LV(i) = 0;
        TMatrix(i,:) = 0;
        TMatrix(i,i) =1 ;
        SPM(i,:) = 0;
  
        
       
    end
end

for j = 1:M+1
    RHS = SPM*Q+LV;
    Q = TMatrix\RHS;
    trisurf(LNODES,x,y,Q(:,:))
%shading interp
xlabel('r','fontsize',16) 
xlim([0 xmax])
ylim([0 xmax])
zlim([0 1])
ylabel('z','fontsize',16)
zlabel('q(r,z)','fontsize',16)
title(['Partial pressure at t= ',num2str(T(j))],'fontsize',20)
pause(1e-10)



   
end




