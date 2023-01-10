function y = AATsolve(rhs,n,Xscale)
   %if Xscale == 1; rhs = rhs/n; end
   rhs = rhs*Xscale^2/n;
   y1 = rhs(1:n) - rhs(1);
   y2 = rhs(n+1:end)-mean(y1);
   y = [y1;y2];
end