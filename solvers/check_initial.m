function yflag = check_initial(y0,G,Ainput);
   ATy0 = Ainput.ATmap(y0);
   tmp = ATy0 + G;
   if max(max(tmp)) < -1e-3
      fprintf('\n initial point is not good');
      fprintf('\n using the intial point generated with out nonegative constraints');
      yflag = 1;
   else
      yflag = 0;
   end
end