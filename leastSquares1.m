%% -*- mode : octave ; -*-
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = leastSquares1 ()
%% function [] = leastSquares1 ()
x = [0: pi /8: pi ];
for k = 0:6
c(k +1) = sum (x.^ k) ;
end
for k = 0:3
b(k +1) = sum (x .^ k .* cos (x)) ;
end
for i = 1:4
for j = 1:4
A(i ,j) = c (i+j -1) ;
end
end
a = A\b';
% polynomial evaluation points .
y = [0:0.01*2* pi :2* pi ];
% evaluate the polynomial at the data points
p = zeros (1 ,101) ;
for i = 0:3
p = p + a (i +1) *y .^ i;
end
hold on
plot (y ,cos (y))
plot (y ,p ,'k')
plot (x ,cos (x) ,'s')