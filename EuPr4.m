clear all;
for i = 1 : 99
    for j = 1 : 99
        if (int2str(i * j) == reverse(int2str(i * j)))
          A(i,j) = i * j;
        else
            A(i,j) = 0;
        end
    end
end

A = A(:)

max(A)