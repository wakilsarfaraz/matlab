clc; clear all
primsum = 0;

for k = 1 : 2000000
    if (isprime(k) == 1)
        primsum = primsum + k;
       
    end
end
primsum