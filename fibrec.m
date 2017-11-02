 function f = fibrec(k)
 if (k==1)
     f= 1;
 elseif (k == 2)
     f = 1;
 else
     f = fibrec(k-1) + fibrec(k-2);   
 end