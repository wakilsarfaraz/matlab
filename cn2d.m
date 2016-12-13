% cn2D.m - Analogue of cn.m but for 2D:
%          the heat equation u_t = u_xx + u_yy on the unit square.
%
%          xx, yy, uu represent data on a 2D mesh.
%          x, y, u represent the same data stretched to 1D vectors,
%             since that's what's needed to do linear algebra.
 
% Grid:
  J = input('J? ');
  dx = 1/J;
  s = (dx:dx:1-dx)';
  dt = .0002;
  nu = dt/dx^2;
  [xx,yy] = meshgrid(s,s);                
  x = xx(:); y = yy(:); 
  J1 = J-1;
 
% Initial condition in the form of a plus sign:
  u = zeros(size(x));
  u(find(abs(x-.5)<.4 & abs(y-.5)<.07)) = 1;
  u(find(abs(x-.5)<.07 & abs(y-.5)<.4)) = 1;
 
% Left-hand-side matrix (could be done more slickly!):
  A = (1+2*nu)*speye(J1^2); 
  for i = 0:J1-1
    for j = 0:J1-2
      A(i*J1+j+1,i*J1+j+2) = -nu/2;
      A(i*J1+j+2,i*J1+j+1) = -nu/2;
      A(j*J1+i+1,j*J1+i+1+J1) = -nu/2;
      A(j*J1+i+1+J1,j*J1+i+1) = -nu/2;
    end
  end
  spy(A), title('sparsity pattern of A')
  pause
   
% Right-hand-side matrix:
  B = (1-2*nu)*speye(J1^2); 
  for i = 0:J1-1
    for j = 0:J1-2
      B(i*J1+j+1,i*J1+j+2) = nu/2;
      B(i*J1+j+2,i*J1+j+1) = nu/2;
      B(j*J1+i+1,j*J1+i+1+J1) = nu/2;
      B(j*J1+i+1+J1,j*J1+i+1) = nu/2;
    end
  end
 
% Time-stepping:
  t = 0;
  while 1
    uu = reshape(u,J1,J1);
    mesh(xx,yy,uu), colormap([0 0 0])  % <- this makes curves black
    axis([0 1 0 1 0 2]), view(-20,70)
    title(['t = ' num2str(t)], 'fontsize',16)
    t = t+dt;
    u = A\(B*u);
    pause
  end