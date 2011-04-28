for i=1:length(maxsz);
   dv = vv(i).meas ./ vr.meas  - 1;
   disp([i, maxsz(i), 1e3*mean(abs(dv)), vv(i).n_ne/1e4]);
end

