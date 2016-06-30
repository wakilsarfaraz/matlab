%% Phase plane for schnakenberg reaction system in absence of diffusion
eps = 0.3;
mul=0.25;
tmax=0.1;
cent = 2;
a = 0.052;
b = 0.276;
m=1;
n=1;
f = @(t,y) [-pi^2*(m^2+n^2)*y(2)+a-y(1)+y(1)^2*y(2);-pi^2*(m^2+n^2)*y(2)+b-y(1)^2*y(2)]
figure(1)
vectfield(f,a+b-mul*eps:.008:a+b+mul*eps,b/(a+b)^2-mul*eps:.008:b/(a+b)^2+mul*eps)
 
  title('Unstable star','fontsize',20)
%   xlabel('u','fontsize',16)
%   ylabel('v','fontsize',16)
hold on

[ts,ys] = ode45(f,[0,tmax],[a+b-mul*eps/cent; b/(a+b)^2]);
 plot(ys(:,1),ys(:,2),'linewidth',2,'color','r')
 
[ts,ys] = ode45(f,[0,tmax],[a+b+mul*eps/cent; b/(a+b)^2]);
 plot(ys(:,1),ys(:,2),'linewidth',2,'color','r')
 
 [ts,ys] = ode45(f,[0,tmax+5],[a+b; b/(a+b)^2+mul*eps/cent]);
 plot(ys(:,1),ys(:,2),'linewidth',2,'color','r')
 
  [ts,ys] = ode45(f,[0,tmax+5],[a+b; b/(a+b)^2-mul*eps/cent]);
 plot(ys(:,1),ys(:,2),'linewidth',2,'color','r')
 
   [ts,ys] = ode45(f,[0,tmax],[a+b+2*mul*eps/(3*cent); b/(a+b)^2+2*mul*eps/(3*cent)]);
 plot(ys(:,1),ys(:,2),'linewidth',2,'color','r')
 
    [ts,ys] = ode45(f,[0,tmax],[a+b+2*mul*eps/(3*cent); b/(a+b)^2-2*mul*eps/(3*cent)]);
 plot(ys(:,1),ys(:,2),'linewidth',2,'color','r')
 
    [ts,ys] = ode45(f,[0,tmax],[a+b-2*mul*eps/(3*cent); b/(a+b)^2+2*mul*eps/(3*cent)]);
 plot(ys(:,1),ys(:,2),'linewidth',2,'color','r')
 
     [ts,ys] = ode45(f,[0,tmax],[a+b-2*mul*eps/(3*cent); b/(a+b)^2-2*mul*eps/(3*cent)]);
 plot(ys(:,1),ys(:,2),'linewidth',2,'color','r')
xlim([(a+b-mul*eps) (a+b+mul*eps)])
ylim([b/(a+b)^2-mul*eps b/(a+b)^2+mul*eps])
xlabel('u','fontsize',18)
ylabel('v','fontsize',18)
 set(findobj('type','legend'),'fontsize',8)
set(findobj('type','axes'),'fontsize',12)
hold off

% figure(2)
% options=odeset('OutputFcn','odephas2');
% ode45(f,[0,200],[1;1],options)
% title('Trajectories in the phase plane')
% xlabel('u')
% ylabel('v')

% 
% figure(3)
% [ts,ys] = ode45(f,[0,20],[1;0]); % find ts, ys, but don't show
% plot(ts,ys) % make plot of y1 and y2 vs. t
% [ts,ys]; % show table with 3 columns for t, y1, y2