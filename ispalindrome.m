function tf = ispalindrome(n)
%this function will check whether the given number is 
%palindrome or not

%Fist we will assume that the given number is palindrome
tf = true;

%Convert the number to string so that we can use the number as
%a vector

n = num2str(n);

%Now we have to find the number of 
%iterations required to check the given number
%If the number is even then length/2
%If the number is odd then (length-1)/2

if rem(length(n),2) == 0
    len = length(n)/2;
else
    len = (length(n)-1)/2;
end

%A for loop to iterate through the number to 
%check for the same numbers

for i = 1:len
    if n(i) ~= n(length(n)+1-i)
        tf = false;
        break;
    end
end