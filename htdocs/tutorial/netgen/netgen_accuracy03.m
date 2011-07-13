for i=1:length(maxsz);
   dv(:,i) = vv(i).meas ./ vr.meas  - 1;
   disp([i, maxsz(i), 1e3*mean(abs(dv(i,:))), vv(i).n_ne/1e4]);
   n_ne(i,:) = vv(i).n_ne;
end

loglog(maxsz,mean(abs(dv)),'*-');
xlabel('Maximum element size');
ylabel('Relative error');
grid 
print_convert netgen_accuracy03a.png '-density 75'

loglog(n_ne(:,1),mean(abs(dv)),'*');
xlabel('Number of nodes');
ylabel('Relative error');
grid on
print_convert netgen_accuracy03b.png '-density 75'

loglog(n_ne(:,2),mean(abs(dv)),'*');
xlabel('Number of elements');
ylabel('Relative error');
grid on
print_convert netgen_accuracy03c.png '-density 75'
