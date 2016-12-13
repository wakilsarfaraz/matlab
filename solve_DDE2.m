% Program to solve a DDE for the protein and mRNA problem with two discrete delays

clear all;
close all;

T=500; dt=0.01;
N=floor(T/dt);

tau_p = 1.7; % time delay in DDE for p
tau_m = 7.1; % time delay in DDE for m

kp = floor(tau_p/dt); % number of steps between -tau and (in general) t
km = floor(tau_m/dt);

% parameter values
a=4.5; 
%a=1.5; % Other choices of 'a' from supplementary material
%a=0.5;
b=0.23; c=0.23;
q=33; pcrit=40; n=2;

% initialising vectors for p and m
p = zeros(kp+1,1);
m = zeros(km+1,1);

% initial conditions for p for each step between -tau_p and 0
for i = 1:kp+1
    p(i) = 0;
end

% initial conditions for m for each step between -tau_m and 0
for i = 1:km+1;
    m(i) = 0;
end

% setting first element of pc and mc for plotting
pc(1) = p(end);
mc(1) = m(end);

for i = 1:N
    pnew = p(end) + (a*m(1)-b*p(end))*dt; % using history to compute next iteration
    mnew = m(end) + (q/(1+(p(1)/pcrit)^n) - c*m(end))*dt;
    
    pn = zeros(kp+1,1);
    mn = zeros(km+1,1);
    
    % pn and mn used to shift p and m along by dt
    for j = 2:kp+1
        pn(j-1) = p(j);
    end
    
    for j = 2:km+1
        mn(j-1) = m(j);
    end
    
    % setting the last element of p and m to pnew and mnew to calculate the
    % next term on the next iteration
    pn(kp+1) = pnew;
    p=pn;
    mn(km+1) = mnew;
    m=mn;
    
    % keeping track of all computed elements of p and m so they can be
    % plotted
    pc(i+1) = p(end);
    mc(i+1) = m(end);
end

time = linspace(0,N*dt,N+1);

plot(time,pc,'b');
hold on
plot(time,mc,'r');
title('Protein molecules and mRNA molecules per cell, with time in minutes');
xlabel('time t');
legend('p(t)','m(t)');

