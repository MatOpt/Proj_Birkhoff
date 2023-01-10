warning off
clear all
addpath(genpath(pwd));
HOME = pwd;
testdatadir = [HOME,filesep,'datafiles'];
eval([' cd ', HOME]); 

fname{1} = 'Proj_2000';



stoptol = 1e-15;
options.stoptol = stoptol;
options.print = 1;

run_Proj_BP = 1;

for kk = 1
   rng('default');
   clear G
   probname = fname{kk}; 
   probtstart = clock;
   load([testdatadir,filesep,probname,'.mat'],'G');
   fprintf('\n load G time = %3.1f', etime(clock, probtstart));
   
   n = length(G);
   
   %construct BP linear constraints
   ONE = ones(n,1);
   b = [ONE;ONE];
   Amap = @(X) [sum(X,2);sum(X,1)'];
   ATmap = @(y) y(1:n)*ONE' + ONE*y(n+1:end)';
   ATAmap = @(X) ATmap(Amap(X));

   Xscale = 1;
   b = b/Xscale;
   Ainput.Amap = @(X) Amap(X)/Xscale;
   Ainput.ATmap = @(y) ATmap(y/Xscale);
   Ainput.ATAmap = @(X) ATAmap(X)/Xscale^2;

   fprintf('\n ===============start proj=====================');
   if run_Proj_BP
      [obj,y,X,info,runhist] = SN_Proj_BP(Ainput,b,n,G,options);
      nal_res = info;
   end
   fprintf('\n ===============end   proj=====================\n');
end
