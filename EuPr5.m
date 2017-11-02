clc; clear all; format long;

smallnumber = 1;

for k = 1 : 20
    if (mod(smallnumber,k+1) ~= 0)
    smallnumber = lcm(smallnumber,(k+1));
    else
        smallnumber = smallnumber;
    end
end
Answer=smallnumber

%Check if it is divisible by all numbers from 1 to 20
remainset = zeros(1, 20);
for i = 1 : 20
    remainset(i) = mod(smallnumber,i);
end

remainset