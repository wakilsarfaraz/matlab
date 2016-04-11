function pdeFEM()

 [p,e,t]=initmesh('squareg');
 [p,e,t]=refinemesh('squareg',p,e,t);

 u0=zeros(size(p,2),1); % initial condition
 ix=find(sqrt(p(1,:).^2+p(2,:).^2)<0.4); % initial condition
 u0(ix)=ones(size(ix)); % initial condition
 tlist=linspace(0,0.1,20); % time interval
 c = 1;
 a = 0;
 f = 0;
 d = 1;
 u1=parabolic(u0,tlist,'squareb1',p,e,t,c,a,f,d);

 figure(1)
 for j=1:10:length(tlist)
 plot3(p(1,:), p(2,:),u1(:,j),'.');
 hold on,
 end

 u0=ones(size(p,2),1); % initial condition
 tlist=linspace(0,1,100); % time interval
 u1=parabolic(u0,tlist,@pdebound,p,e,t,1,0,0,1);

figure(2)
plot3(p(1,:), p(2,:),u1(:,length(tlist)),'.');

 u0=ones(size(p,2),1); % initial condition
 kdeg =2; % degradation rate
 tlist=linspace(0,1,100); % time interval
 u1=parabolic(u0,tlist,@pdebound,p,e,t,1,kdeg,0,1);

 figure(3)
 plot3(p(1,:), p(2,:),u1(:,length(tlist)),'.');
 end

 function [qmatrix,gmatrix,hmatrix,rmatrix] = pdebound(p,e,u,time)

 ne = size(e,2); % number of edges
 qmatrix = zeros(1,ne);
 gmatrix = qmatrix;
 hmatrix = zeros(1,2*ne);
 rmatrix = hmatrix;

 for k = 1:ne
 x1 = p(1,e(1,k)); % x at first point in segment
 x2 = p(1,e(2,k)); % x at second point in segment
 xm = (x1 + x2)/2; % x at segment midpoint
 y1 = p(2,e(1,k)); % y at first point in segment
 y2 = p(2,e(2,k)); % y at second point in segment
 ym = (y1 + y2)/2; % y at segment midpoint

 switch e(5,k)
 case {1} % pick one boundary
 hmatrix(k) = 1;
 hmatrix(k+ne) = 1;
 rmatrix(k) = 1;
 rmatrix(k+ne) = 1;
 otherwise % other boundaries
 qmatrix(k) = 0;
 gmatrix(k) = 0;
 end
 end
 end
