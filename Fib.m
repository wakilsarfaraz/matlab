%% This script returns the Nth element of fibonacci sequence for a given N.
% Name
% Date

clc; clear all;
tic
N = 15 ;

fib = zeros(N,1);

fib(1) = 1;
fib(2) = 1;

for k = 2 : N-1
    fib(k+1) = fib(k)+fib(k-1);
end
fib
max(fib)
toc
