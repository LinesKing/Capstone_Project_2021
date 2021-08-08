function [ac] = accept(x,xa,xL,xU)
% parameters
filterOn = false;
betterOn = false;
n = size(x,2);
gammaTh = 10^-5;
gammaPhi = 10^-5;

% check for violations of variable constraints
% shouldn't need to check these
viol = 0;
for i = 1:n
    if(xa(i)<xL(i))
        viol = viol + 1;
    end
    if(xa(i)>xU(i))
        viol = viol + 1;
    end
end

if (viol==0)
  % either condition better compared to last iteration
  if (betterOn)
     better = (phi(bp,xa,xL,xU) <= phi(bp,x,xL,xU) - gammaPhi*bp_theta(bp,x)) || (theta(bp,xa) <= (1-gammaTh)* theta(bp,x));
  else
     better = true;
  end
     
  % apply filter to determine whether to accept
  if (better)
      ac = true;
      if (filterOn)
         nf = size(filter,1);
         for i = 1:nf
             if (((1-gammaTh)*bp_theta(bp,xa) > filter(i,1)) || (bp_phi(bp,xa,xL,xU)-gammaPhi*bp_theta(bp,xa) > filter(i,2))),
                 ac = false;
             end
         end
      end
  else
      ac = false;
  end
else
    ac = false;
end
