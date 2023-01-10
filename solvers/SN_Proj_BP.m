%%*************************************************************************
%% Using semimsooth Newton method to solve the following
%% inf 0.5*||Pi_{>=0}(Aty + G)||^2 - 0.5*||G||^2 - b'*y
%%
%% modifyied by LXD on Jan 11, 2017 to include initial scaling
%%*************************************************************************
function [obj,y,X,info,runhist] = SN_Proj_BP(Ainput,b,n,G,options,y0)


print = 1;
breakyes = 0;
maxiter = 500;
tiny = 1e-10;
stoptol = 1e-9;
maxitpsqmr =500;
precond = 1;
scale = 1;
initial_op = 2;
alpha = 1;
alpha_scale = 0;
infomore = 0;
runhist.psqmr = 0;
runprojalone = 0;

if isfield(options,'print'); print = options.print; end
if isfield(options,'maxiter'); maxiter = options.maxiter; end
if isfield(options,'tiny'); tiny = options.tiny; end
if isfield(options,'stoptol'); stoptol = options.stoptol; end
if isfield(options,'precond'); precond = options.precond; end
if isfield(options,'scale'); scale = options.scale;end
if isfield(options,'infomore'); infomore = options.infomore;end
if isfield(options,'runprojalone'); runprojalone = options.runprojalone; end


fprintf('\n*************************************************************************************');
fprintf('\n*  SN_Proj_BP: A dual semismooth Newton method for the Birkhoff polytope projection *');
fprintf('\n*  Authors: Xudong Li, Defeng Sun, and Kim-Chuan Toh                                *');
fprintf('\n*************************************************************************************');

%% preperation    
tstart = clock;
m = length(b);
borg = b;
normborg = 1 + norm(borg);
%normb = normborg;
Gorg = G;
FnormGorg = mexFnorm(Gorg);
%normGorg = 1 + FnormGorg;
FnormG = FnormGorg;
par.n = n;
if norm(G - speye(n,n),'fro') < stoptol
   fprintf('\n at this level, I is a good solution');
   obj = 0;
   X = speye(n,n);
   y = zeros(m,1);
   info.AX = ones(2*n,1);
   info.ATy = zeros(n,n);
   info.P = speye(n,n);
   info.Pe = ones(n,1);
   info.eTP = ones(n,1);
   runhist = [];
   return
end
if isstruct(Ainput)
   if isfield(Ainput,'Amap'); Amap0 = Ainput.Amap; end
   if isfield(Ainput,'ATmap'); ATmap0 = Ainput.ATmap; end
end
if ~exist('y0','var')
   if max(max(G)) < -1e-3 && initial_op == 1; initial_op = 2;end
   %initial_op = 1;
   if initial_op == 1
      y = zeros(m,1);  
   elseif initial_op == 2
      rhs_initial = b - Amap0(G);
      y = AATsolve(rhs_initial,n,1);
   end
else
   y_regenerate = check_initial(y0,G,Ainput);
   if y_regenerate
      rhs_initial = b - Amap0(G);
      y = AATsolve(rhs_initial,n,1);
   else
      y = y0; 
   end
end


if scale == 1
   Xscale = sqrt(n); % max(max(abs(G)));%max(max(abs(G))); %max(max(abs(G))); %sqrt(n);
   b = b/Xscale;
   normb = 1 + norm(b);
   Amap = @(x) Amap0(x)/Xscale;
   ATmap = @(y) ATmap0(y/Xscale);
   y = Xscale*y;
else
   Xscale = 1;
   Amap = Amap0;
   ATmap = ATmap0;
end

maxG = max(max(abs(G)));
if maxG > 2; alpha_scale = 1; end
if alpha_scale == 1
   alpha = max(1,maxG);
   Amap = @(X) Amap(X)*alpha;
   ATmap = @(y) ATmap(y*alpha);
   G = G/alpha;
   FnormG = FnormG/alpha;
   y = y/alpha^2;
end
par.scale = scale;
par.Xscale = Xscale;
par.alpha = alpha;
Ainput_nal.Amap = Amap;
Ainput_nal.ATmap = ATmap;
matvecfname = 'matvecAPAT1118';%'matvecAPAT'; 
%%
Aty = ATmap(y);
Sinput = Aty + G;
X = max(Sinput,0);
P = (Sinput > 0); 
spP = sum(sum(P))/n^2;
par.nP = [];
if (spP < 0.3)
   P = sparse(P);
   par.spP = 1;
elseif (spP > 0.7)
   par.spP = -1;
   par.nP = sparse(~P);
else
   par.spP = 0;
end  
par.P = P;
par.Pe = sum(P,2);
par.eTP= sum(P,1)';
AX = Amap(X);
Rp = AX - b;
step_op.FnormG = FnormG;
Ly = b'*y - 0.5*mexFnorm(X)^2 +0.5*FnormG^2;
%runhist.psqmr(1) = 0;
%runhist.findstep(1) = 0;
if print
   fprintf('\n Using SSNCG to solve the (D) problem , with initial_op = %d', initial_op);
   fprintf('\n problem size n = %3.0f, m = %3.0f', n, m);
   fprintf('\n problem scale Xscale = %3.2e, alpha = %3.2e',Xscale,alpha);
   fprintf('\n ----------------------------------------------------------');
   fprintf('\n iter      dobj       dgrad    time  | cg:tol  res      NO.');
   fprintf('     | l_step    l_it');
end
epsilonop = 0;
findstep_op = 1;
iterstep = 0;
%% mian Newton iteration
for itersub = 1:maxiter    
%    yold = y; Atyold = Aty;
    GradLy = -Rp;
    normGradLy = norm(GradLy)*Xscale/normborg;
    priminf_sub = normGradLy; 
    runhist.priminf(itersub) = priminf_sub;
    runhist.Ly(itersub)      = Ly;
    ttime = etime(clock,tstart);
    if (print)
        fprintf('\n%2.0d  %- 11.10e %3.2e  %3.1f',...
                 itersub,alpha^2*Ly,priminf_sub,ttime);
    end
    if max([normGradLy]) < max(stoptol) 
        msg = 'Problem solved';
        if print
            fprintf('\n%s  ',msg);
            fprintf(' gradLy = %3.2e, stoptol=%3.2e, time = %3.1f',normGradLy,stoptol,ttime);
        end
        breakyes = -1;
        break;
    end
    %% Compute Newton direction
     if epsilonop == 1
        if normGradLy > 1
           par.epsilon = 0.01*normGradLy;
        elseif normGradLy > 1e-2
           par.epsilon = min(1e-2,0.1*normGradLy);
        elseif normGradLy > 1e-4
           par.epsilon = min(1e-4,0.1*normGradLy);
        end
     elseif epsilonop == 0
        par.epsilon = min(1e-4,0.1*normGradLy);
     elseif epsilonop == 2
        par.epsilon = min(1e-2,0.01*normGradLy);
     end
      %% good to add
     par.precond = precond;
     if precond == 1
        par.invdiagM = [1./(par.Pe + par.epsilon);1./(par.eTP+par.epsilon)];
     end
     if normGradLy > 1
        maxitpsqmr = 200;
     else
        maxitpsqmr = 200;
     end
     if (itersub > 1) 
          prim_ratio = priminf_sub/runhist.priminf(itersub-1); 
     else
          prim_ratio = 0; 
     end
     rhs = GradLy;
     normrhs = norm(rhs);
     tolpsqmr = min([5e-3,1e-4*normrhs]);
%      if normrhs > 1
%         tolpsqmr = min(5e-3, 0.05*normrhs);
%      elseif normrhs > 1e-2
%         tolpsqmr = min(1e-3, 0.05*normrhs);
%      elseif normrhs > 1e-4
%         tolpsqmr = min(5e-4, 0.01*normrhs);
%      else
%         tolpsqmr = min(1e-4, 0.1*normrhs);
%      end
     const2 = 1;
     if itersub > 1 && (prim_ratio > 0.5 || priminf_sub > 0.1*runhist.priminf(1))
        const2 = 0.5*const2;
     end
     tolpsqmr = const2*tolpsqmr;
     par.tol = tolpsqmr; par.maxit = maxitpsqmr;
     [dy,~,resnrm,solve_ok] =  psqmry(matvecfname,Ainput_nal,rhs,par); 
     Atdy = ATmap(dy);
     iterpsqmr = length(resnrm)-1;
     if iterpsqmr ==0; keyboard; end
     if (print)
          fprintf('  | %3.1e %3.1e %3.0d',par.tol,resnrm(end),iterpsqmr);
          fprintf(' %2.1f',const2);
     end
     par.iter = itersub;
     if (itersub <=3) 
         stepop = 1;
     else
         stepop = 2;
     end
%      if priminf_sub > 1e-2; 
%         step_op.alphamax = 1.99;
%      else
%         step_op.alphamax = 1;
%      end
     steptol = 1e-5; step_op.stepop=stepop;
     if iterstep == 17; findstep_op = 0; end %% important: LXD
     if findstep_op == 1
        [par,Ly,y,Aty,X,AX,alp,iterstep] = ...
           findstep(Ainput_nal,par,b,G,Ly,y,Aty,X,AX,dy,Atdy,steptol,step_op); 
     else
        [par,Ly,y,Aty,X,alp,iterstep] = ...
           findstep_old(par,b,G,Ly,y,Aty,X,dy,Atdy,steptol,step_op);
        AX = Amap(X);
     end
     %AX = Amap(X);
     Rp = AX- b;
     runhist.solve_ok(itersub) = solve_ok;
     runhist.psqmr(itersub)    = iterpsqmr; 
     runhist.findstep(itersub) = iterstep; 
     %if alp < tiny; breakyes =11; break; end
     Ly_ratio = 1; 
     if (itersub > 1)
          Ly_ratio = (Ly-runhist.Ly(itersub-1))/(abs(Ly)+eps);
     end
     if (print)
          fprintf(' | %3.2e %- 2.0f',alp,iterstep);
          if (Ly_ratio < 0) && (-Ly_ratio > 1e-6); fprintf('-'); end
     end
end

info.time  = etime(clock,tstart);
info.maxCG = max(runhist.psqmr);
info.totalCG = sum(runhist.psqmr);
info.avgCG = info.totalCG/itersub;
info.breakyes = breakyes;
info.itersub = itersub;
info.eta = priminf_sub;
info.iter= itersub;
y = y/Xscale;

obj(1) = 0.5*mexFnorm(X - G)^2;
obj(2) = Ly;
if alpha ~= 1
   y = y*alpha^2;
   Aty = alpha*Aty;
   obj(2) = alpha^2*obj(2);
   X = alpha*X;
   obj(1) = alpha^2*obj(1);
   if print
      fprintf('\n  primal obj_val = %11.10e',obj(1));
      fprintf('\n    dual obj_val = %11.10e',obj(2)); 
   end
end
info.obj = obj;
gap = obj(1) - obj(2); 
rel_gap = gap/(1 + abs(obj(1)) + abs(obj(2)));

info.gap = gap;
info.rel_gap = rel_gap;
info.eta_hist = runhist.priminf;
info.n = n;

if runprojalone
   normX = mexFnorm(X);
   info.etaXorg = mexFnorm(max(-X,0));
   info.etaX = info.etaXorg/(1 + normX);
   info.etaAorg = norm(Amap0(X) - borg);
   info.etaA = info.etaAorg/normborg;
   info.etaPorg = max(info.etaAorg, info.etaXorg);
   info.etaP = max(info.etaA,info.etaX);
end

if infomore
   info.AX = AX*Xscale;
   info.ATy = Aty;
   info.P = par.P;
   info.Pe = par.Pe;
   info.eTP = par.eTP;
   info.spP = par.spP;
   info.nP  = par.nP;
end

if print
   fprintf('\n gap = %3.2e, rel_gap = %3.2e', gap, rel_gap); 
   if runprojalone
      fprintf('\n etaXorg = %3.2e,  etaX = %3.2e ', info.etaXorg, info.etaX);
      fprintf('\n etaAorg = %3.2e,  etaA = %3.2e ', info.etaAorg, info.etaA);
      fprintf('\n etaPorg = %3.2e,  etaP = %3.2e ', info.etaPorg, info.etaP);
   end
end
if print
   fprintf('\n Newton step total CG = %3.0d, CG per iter = %3.1f',info.totalCG, info.totalCG/itersub);
   fprintf('\n ---------------End Newton method ----------------------\n');
end
end
%%********************************************************************
function [par,Ly,y,Aty,X,AX,alp,iter] = ...
         findstep(Ainput,par,b,G,Ly0,y0,Aty0,X0,AX0,dy,Atdy,tol,options)
   alphamax = 1;
   if isfield(options,'FnormG'); FnormG = options.FnormG; end
   if ~exist('FnormG','var'); FnormG = norm(G,'fro'); end
   if isfield(options,'alphamax'); alphamax = options.alphamax; end
   printlevel = 0; 
   maxit = ceil(log(1/(tol+eps))/log(2));
   c1 = 1e-4; c2 = 0.9; 
   n = par.n;
%%
   dytb = dy'*b;
   y0tb = y0'*b;
   g0  = dytb - dy'*AX0; %sum(sum(Atdy.*X0));%mexSum2AB(Atdy,sparse(X0));%
   if (g0 <= 0)
      alp = 0; iter = 0; 
      if (printlevel) 
         fprintf('\n Need an ascent direction, %2.1e  ',g0); 
      end
      y = y0;
      Aty = Aty0;
      AX = AX0;
      X = X0;
      Ly = Ly0;
      return;
   end  
   breakyes = 0;
%%
   %alp = 1; 
   alpconst = 0.5;
   Aty0pG = Aty0 + G;
   for iter = 1:maxit
      if (iter==1)          
         alp = alphamax; 
      else
         alp = alpconst*alp;
      end
      y = y0 + alp*dy;
      %Aty = Aty0 + alp*Atdy;
      Sinput = alp*Atdy + Aty0pG;
      X = max(Sinput,0);
      Ly   = y0tb + alp*dytb - 0.5*mexFnorm(X)^2 + 0.5*FnormG^2;
      if printlevel
          fprintf('\n ------------------------------------- \n');
          fprintf('\n alp = %2.2f, LQ = %11.10e, LQ0 = %11.10e',alp,Ly,Ly0);
          fprintf('\n ------------------------------------- \n');
      end
      if (iter==1)
         AX = Ainput.Amap(X);
         galp = dytb - dy'*AX; %sum(sum(Atdy.*X));
         gLB = g0; gUB = galp; 
         if (sign(gLB)*sign(gUB) > 0)
            if (printlevel); fprintf('|'); end
            Aty = Aty0 + alp*Atdy;
            P = Sinput > 0;
            spP = sum(sum(P))/n^2;
            par.nP = [];
            if (spP < 0.3)
               P = sparse(P);
               par.spP = 1;
            elseif (spP > 0.7)
              par.spP = -1;
              par.nP = sparse(~P);
            else
              par.spP = 0;
            end  
            %P = sparse(Sinput > 0);
            par.P = P;
            par.Pe = sum(P,2);
            par.eTP= sum(P,1)';
            %return;             
            breakyes = 1;
         end
      end
      if (Ly-Ly0-c1*alp*g0 > -1e-8/max(1,abs(Ly0))) && breakyes == 0
            Aty = Aty0 + alp*Atdy;
            P = Sinput > 0;
            spP = sum(sum(P))/n^2;
            par.nP = [];
            if (spP < 0.3)
               P = sparse(P);
               par.spP = 1;
            elseif (spP > 0.7)
              par.spP = -1;
              par.nP = sparse(~P);
            else
              par.spP = 0;
            end
            par.P = P;
            par.Pe = sum(P,2);
            par.eTP= sum(P,1)';
            if iter >1; AX = Ainput.Amap(X); end
            %return;
            breakyes = 2;
      end
      if (alp < 1) && printlevel
          fprintf('\n iter = %2d, ------line search value------------\n',iter);
          fprintf('\n ------alp = %2.2f, LQ = %11.10e, LQ0 = %11.10e',alp,Ly,Ly0);       
      end
      if breakyes > 0; break; end
   end 
   if iter == maxit && breakyes == 0
      Aty = Aty0 + alp*Atdy;
      P = Sinput > 0;
      spP = sum(sum(P))/n^2;
      par.nP = [];
      if (spP < 0.3)
         P = sparse(P);
         par.spP = 1;
      elseif (spP > 0.7)
         par.spP = -1;
         par.nP = sparse(~P);
      else
         par.spP = 0;
      end
      par.P = P;
      par.Pe = sum(P,2);
      par.eTP= sum(P,1)';
      AX = Ainput.Amap(X);
   end
   if (printlevel); fprintf('m'); end
end
%%********************************************************************

function [par,Ly,y,Aty,X,alp,iter] = ...
         findstep_old(par,b,G,Ly0,y0,Aty0,X0,dy,Atdy,tol,options)
   if isfield(options,'stepop'); stepop = options.stepop; end
   if isfield(options,'FnormG'); FnormG = options.FnormG; end
   if ~exist('FnormG','var'); FnormG = norm(G,'fro'); end
   printlevel = 0; 
   maxit = ceil(log(1/(tol+eps))/log(2));
   c1 = 1e-4; c2 = 0.9; 
   n = par.n;
%%
   dytb = dy'*b;
   y0tb = y0'*b;
   g0  = dytb - sum(sum(Atdy.*X0));
   if (g0 <= 0)
      alp = 0; iter = 0; 
      if (printlevel) 
         fprintf('\n Need an ascent direction, %2.1e  ',g0); 
      end
      y = y0;
      Aty = Aty0;
      X = X0;
      Ly = Ly0;
      return;
   end  
%%
   alp = 1; alpconst = 0.5;
   Aty0pG = Aty0 + G;
   for iter = 1:maxit
      if (iter==1)          
         alp = 1; LB = 0; UB = 1; 
      else
         alp = alpconst*(LB+UB);
      end
      y = y0 + alp*dy;
      %Aty = Aty0 + alp*Atdy;
      Sinput = alp*Atdy + Aty0pG;
      X = max(Sinput,0);
      galp = dytb - sum(sum(Atdy.*X));
      Ly   = y0tb + alp*dytb - 0.5*mexFnorm(X)^2 + 0.5*FnormG^2;
      if printlevel
          fprintf('\n ------------------------------------- \n');
          fprintf('\n alp = %2.2f, LQ = %11.10e, LQ0 = %11.10e',alp,Ly,Ly0);
          fprintf('\n ------------------------------------- \n');
      end
      if (iter==1)
         gLB = g0; gUB = galp; 
         if (sign(gLB)*sign(gUB) > 0)
            if (printlevel); fprintf('|'); end
            Aty = Aty0 + alp*Atdy;
            P = sparse(Sinput > 0);
            spP = sum(sum(P))/n^2;
            par.nP = [];
            if (spP < 0.3)
               P = sparse(P);
               par.spP = 1;
            elseif (spP > 0.7)
               par.spP = -1;
               par.nP = sparse(~P);
            else
               par.spP = 0;
            end
            par.P = P;
            par.Pe = sum(P,2);
            par.eTP= sum(P,1)';
            return;             
         end
      end
      if (abs(galp) < c2*abs(g0)) && (Ly-Ly0-c1*alp*g0 > -1e-8/max(1,abs(Ly0)))
         if (stepop==1) || ((stepop == 2) && (abs(galp) < tol))
            if (printlevel); fprintf(':'); end
            Aty = Aty0 + alp*Atdy;
            P = sparse(Sinput > 0);
            spP = sum(sum(P))/n^2;
            par.nP = [];
            if (spP < 0.3)
               P = sparse(P);
               par.spP = 1;
            elseif (spP > 0.7)
               par.spP = -1;
               par.nP = sparse(~P);
            else
               par.spP = 0;
            end
            par.P = P;
            par.Pe = sum(P,2);
            par.eTP= sum(P,1)';
            return;
         end
      end
      if (sign(galp)*sign(gUB) < 0)
         LB = alp; gLB = galp;
      elseif (sign(galp)*sign(gLB) < 0) 
         UB = alp; gUB = galp; 
      end
      if (alp < 1) && printlevel 
          fprintf('\n iter = %2d, ------line search value------------\n',iter);
          fprintf('\n ------alp = %2.2f, LQ = %11.10e, LQ0 = %11.10e',alp,Ly,Ly0);       
      end
   end 
   if iter == maxit
      Aty = Aty0 + alp*Atdy;
      P = sparse(Sinput > 0);
      spP = sum(sum(P))/n^2;
      par.nP = [];
      if (spP < 0.3)
         P = sparse(P);
         par.spP = 1;
      elseif (spP > 0.7)
         par.spP = -1;
         par.nP = sparse(~P);
      else
         par.spP = 0;
      end
      par.P = P;
      par.Pe = sum(P,2);
      par.eTP= sum(P,1)';
   end
   if (printlevel); fprintf('m'); end
%%********************************************************************
end
% 
