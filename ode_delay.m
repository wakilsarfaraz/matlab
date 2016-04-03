clear
% parameters
epsilon = 1; % amplitude
kappa = 0; % damping
stx = 200; % steps of 
sty = 100; % steps of b
b_st = -1.5; % starting value for b
b_fi = 1.5; % final value for b
delta_st = -1; % starting value for 
delta_fi = 5; % final value for 
% computational parameters
k = 40; % number of discretization intervals over one period
 intk = 20; % number of integration interval in each step
T = 2*pi; % time period
dt = T/k; % discretization interval length
ddt = dt/intk; % integration interval length
tau = 2*pi; % time delay
m = floor((tau + dt/2)/dt);
wa = (m*dt + dt/2 -tau)/dt;
wb = (tau + dt/2 - m*dt)/dt;
D = zeros(m + 2, m + 2); % matrix D
d = ones(m + 1, 1);
d(1 : 2) = 0;
D = D + diag(d, -1);
D(3, 1) = 1;
for i = 1 : k
c(i) = sum(cos(((i - 1)*2*pi/k) : ddt : ((i -1)*2*pi/k + ddt*(intk - 1))))/intk;
end % discrete values of ci
% start of computation
for y = 1 : sty + 1 % loop for b
b = b_st + (y - 1)*(b_fi -b_st)/sty;
for x = 1 : stx + 1 % loop for 
delta = delta_st + (x - 1)*(delta_fi - delta_st)/stx;
Fi = eye(m + 2, m + 2);
% construct transition matrix Fi
for i = 1 : k
A = zeros(2, 2); % matrix Ai
A(1, 2) = 1;
A(2, 1) = -delta - epsilon*c(i);
A(2, 2) = -kappa;
B = zeros(2, 2); % matrix B
B(2, 1) = b;
P = expm(A*dt); % matrix Pi
R = (expm(A*dt) -eye(2))*inv(A)*B; % matrix Ri
D(1 : 2, 1 : 2) = P;
D(1 : 2, m + 1) = wa*R(1 : 2, 1 : 1);
D(1 : 2, m + 2) = wb*R(1 : 2, 1 : 1);
Fi = D*Fi; % Floquet transition matrix 
end
delta_m(x, y) = delta; % matrix of 
b_m(x, y) = b; % matrix of b
ei_m(x, y) = max(abs(eig(Fi))); % matrix of eigenvalues
end
sty + 1 -y % counter
end
figure
contour(delta_m,b_m,ei_m,[1, 1],'k')