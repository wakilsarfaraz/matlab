%% Phase plane for odes
eps = 0.1;
d=0.3;
tmax=1;
cent = 2;
a = 0.1;
b = 2;
f = @(t,y) [a-y(1)+y(1)^2*y(2);b-y(1)^2*y(2)]
figure(1)
vectfield(f,a+b-d*eps:.008:a+b+d*eps,b/(a+b)^2-d*eps:.008:b/(a+b)^2+d*eps)
 
  title('Unstable star','fontsize',20)
%   xlabel('u','fontsize',16)
%   ylabel('v','fontsize',16)
hold on

[ts,ys] = ode45(f,[0,tmax],[a+b-d*eps/cent; b/(a+b)^2]);
 plot(ys(:,1),ys(:,2),'linewidth',2,'color','r')
 
[ts,ys] = ode45(f,[0,tmax],[a+b+d*eps/cent; b/(a+b)^2]);
 plot(ys(:,1),ys(:,2),'linewidth',2,'color','r')
 
 [ts,ys] = ode45(f,[0,tmax+5],[a+b; b/(a+b)^2+d*eps/cent]);
 plot(ys(:,1),ys(:,2),'linewidth',2,'color','r')
 
  [ts,ys] = ode45(f,[0,tmax+5],[a+b; b/(a+b)^2-d*eps/cent]);
 plot(ys(:,1),ys(:,2),'linewidth',2,'color','r')
 
   [ts,ys] = ode45(f,[0,tmax],[a+b+2*d*eps/(3*cent); b/(a+b)^2+2*d*eps/(3*cent)]);
 plot(ys(:,1),ys(:,2),'linewidth',2,'color','r')
 
    [ts,ys] = ode45(f,[0,tmax],[a+b+2*d*eps/(3*cent); b/(a+b)^2-2*d*eps/(3*cent)]);
 plot(ys(:,1),ys(:,2),'linewidth',2,'color','r')
 
    [ts,ys] = ode45(f,[0,tmax],[a+b-2*d*eps/(3*cent); b/(a+b)^2+2*d*eps/(3*cent)]);
 plot(ys(:,1),ys(:,2),'linewidth',2,'color','r')
 
     [ts,ys] = ode45(f,[0,tmax],[a+b-2*d*eps/(3*cent); b/(a+b)^2-2*d*eps/(3*cent)]);
 plot(ys(:,1),ys(:,2),'linewidth',2,'color','r')
xlim([(a+b-d*eps) (a+b+d*eps)])
ylim([b/(a+b)^2-d*eps b/(a+b)^2+d*eps])
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