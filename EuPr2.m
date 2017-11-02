clear all;
for k = 1 : largefib(4000000)
    A(k) = largefib(k);
    B(k) = largefib(k);
end
C = union(A,B);

for k = 1 : length(C)
    if (mod(C(k),2)==0) 
        C(k) = C(k);
    else
        C(k) = 0;
    end
end
Answer = sum(C)