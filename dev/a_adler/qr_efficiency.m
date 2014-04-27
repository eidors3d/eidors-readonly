N_elec = 16; N_tries = 10;
img = mk_image( ng_mk_cyl_models([2,1,.1],[N_elec,1],[0.1]), 1 );
%for N_pat =  [2,5,10,20,50,100,200,500,1000,2000,5000,10000]
for N_pat =  [500, 2000]
   time = 0;
   for i=1:N_tries % Tries
      patterns = ceil( N_elec*rand(N_pat,4) );
      img.fwd_model.stimulation = stim_meas_list(patterns,16);
      t= cputime;
        vh= fwd_solve(img);
      time = time + (cputime - t);
   end
   fprintf('Patterns = %5d  time = %4.2f \n', [N_pat, time/N_tries]);
end
