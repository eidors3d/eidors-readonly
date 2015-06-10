d = dir;
idx = cell2mat({d(:).isdir});
d = d(idx);

for i = 1:numel(d)
   if strcmp(d(i).name,'..')
      continue
   end
%    fprintf('%s\n', d(i).name);
   cd(d(i).name);
   p = which('eidors_var_id');
   try
      id = eidors_var_id([]);
      if strcmp(id,'id_DA39A3EE5E6B4B0D3255BFEF95601890AFD80709');
         fprintf('%s: CORRECT\n', p);
      else
         fprintf('%s: WRONG\n', p);
      end
   catch
      fprintf('%s: ERROR\n', p);
   end
   if ~strcmp(d(i).name,'.')
      cd ..
   end
end