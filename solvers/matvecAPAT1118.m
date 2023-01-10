%%**************************************
%% matrix vector multiplication 
%% using the sparsity of P
%%**************************************
function y = matvecAPAT1118(x,par,Ainput)
x1 = x(1:par.n);
x2 = x(par.n+1:end);
if par.spP >= 0
   y1 = par.Pe.*x1 + par.P*x2;
   y2 = (x1'*par.P)' + par.eTP.*x2;
elseif par.spP == -1
   y1 = par.Pe.*x1 + (sum(x2)*ones(length(x2),1) - par.nP*x2);
   y2 = (sum(x1)*ones(1,length(x1)) - x1'*par.nP)' + par.eTP.*x2;
end
y = [y1;y2];
if par.scale == 1
   y = y/par.Xscale^2;
end
y = par.alpha^2*y;
%y = Ainput.Amap(par.P.*Ainput.ATmap(x));
if isfield(par,'epsilon')
   y = y + par.epsilon*x;
end
