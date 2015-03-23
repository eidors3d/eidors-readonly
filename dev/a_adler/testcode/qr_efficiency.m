N_elec = 16; N_tries = 2;
%img = mk_image( ng_mk_cyl_models([2,1,.1],[N_elec,1],[0.1]), 1 );
 img = mk_image( mk_common_model('k2c2',N_elec) );
%for N_pat =  [2,5,10,20,50,100,200,500,1000,2000,5000,10000]
%for N_pat =  [50, 500, 5000]
 for N_pat =  500;
   time = 0;
   for i=1:N_tries % Tries
      patterns = ceil( N_elec*rand(N_pat,4) );
      for j=1:N_pat;
         if diff(patterns(j,1:2))==0; patterns(j,1:2)=[1,N_elec]; end
         if diff(patterns(j,3:4))==0; patterns(j,3:4)=[1,N_elec]; end
      end
      img.fwd_model.stimulation = stim_meas_list(patterns,16);
      t= cputime;
      if 0
        vh= fwd_solve(img);
      else
        J= calc_jacobian(img);
      end
      time = time + (cputime - t);
   end
   fprintf('Patterns = %5d  time = %4.2f \n', [N_pat, time/N_tries]);
end
testdata=[
   10   0.57 
   20   0.64 
   50   0.99 
  100   1.44 
  200   1.17 
  500   1.51 
 1000   2.29 
 2000   4.40 
 5000   7.77 
10000   16.55 ];

