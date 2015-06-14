SZ = 1024^3;
eidors_cache('cache_size',SZ);


eidors_msg log_level 0; % silence

   fid = fopen('cache_test.txt','w');
   fprintf(fid,'size\t1st\t2nd\n');
   fclose(fid);

for i = [1:5 10 20 50 100 200 500 1000 2000 5000 10000]
   eidors_cache clear
   eidors_cache('cache_size',round(SZ/i));
   fprintf('CACHE_SIZE:\t %d b\n', round(SZ/i));
     
   t0 = tic;
   inv_solve_abs_core UNIT_TEST
   T1 = toc(t0);
   close all
  
   t0 = tic;
   inv_solve_abs_core UNIT_TEST
   T2 = toc(t0);
   close all
  
   fprintf('1st iteration:\t %f s\n',T1)
   fprintf('2nd iteration:\t %f s\n',T2)
   
   
   fid = fopen('cache_test.txt','a');
   fprintf(fid,'%d\t%f\t%f\n', ...
      eidors_cache('cache_size'), ...
      T1,T2);
   fclose(fid);
end