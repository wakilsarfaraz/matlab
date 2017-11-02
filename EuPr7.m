
clc; clear all

nprim = 1;

nthnum = 10001;

while (length(primes(nprim)) < nthnum)
    nprim = nprim +1;
end

Answer = nprim

