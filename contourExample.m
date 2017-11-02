%Example:
          x=linspace(-3,3,39);
          y=linspace(-2.5,2.5,49);
          [xx,yy]=meshgrid(x,y);
          zz=peaks(xx,yy);
          v=-3:1:5; % contour levels
          subplot(1,2,1)
          [C,h]=contour(xx,yy,zz,v);   % standard contour for comparison
          %clabel(C)
          title Contour

          idx=randperm(numel(zz));     % grab some scattered indices
          n=idx(1:ceil(numel(zz)/2))'; % one half of them
          x=xx(n);                     % get scattered data
          y=yy(n);
          z=zz(n);
          tri=delaunay(x,y);           % triangulate scattered data
          subplot(1,2,2)
          [C,h]=tricontour(tri,x,y,z,v);
          %clabel(C,h)
          title TriContour