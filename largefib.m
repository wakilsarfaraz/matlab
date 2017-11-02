function lfib = largefib(x)
% compute the largest fibonacci number <= x
% Code only works for positive values of x, even though 
% the Fibonacci sequence can be extended for a
% negative index.
if x <= 1
  lfib = 1;
  return
end
fnminus2 = 0;
fnminus1 = 1;
fn = -inf;
while fn <= x
  fn = fnminus1 + fnminus2;
  fnminus2 = fnminus1;
  fnminus1 = fn;
end
lfib = fnminus2;