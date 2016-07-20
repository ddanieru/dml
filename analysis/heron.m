#!/usr/bin/octave -q

function x = heron(x,y)
   x = .5*(x + y/x);
end

function x = errheron(x,y)
   printf("Value     Error      Rel. Error\n");
   printf("-------  -------    ------------\n");
   for i=1:5
      x = heron(x,y);
      err = x - sqrt(y);
      relerr = (x/sqrt(y)) - 1.;
      printf("%7.4f  %7.4e  %7.4e \n",x, err,relerr);
   end
end

# Square root of 2 starting in 1.
x = 1.
y = 2.

errheron(x,y);
