function eidors_deprecate(oldname, newname)

if ischar(oldname) && strcmp(oldname,'UNIT_TEST'), do_unit_test, return, end

eidors_dir = fileparts(which('eidors_obj'));
oldfile = which(oldname, '-all');
if length(oldfile) > 1
   id = find(cellfun('isempty', strfind(oldfile,'deprecated')));
   oldfile = oldfile(id(1));
end
olddir     = fileparts(oldfile{1});
if ~strcmp(olddir(1:min(length(olddir),length(eidors_dir))),eidors_dir)
   eidors_msg('Cannot locate %s in EIDORS',oldname,1);
   return
end

oldpath = sprintf('%s/%s.m',olddir, oldname);
newpath = sprintf('%s/deprecated/%s.m',eidors_dir,oldname);




if nargin == 1
   copy_warn(oldpath, newpath);
   system(sprintf('rm %s',oldpath));
else
   copy_warn(oldpath, newpath, newname);
   % rename all references 
   eidors_rename(oldname, newname);
end



function copy_warn(oldpath, newpath, newname)
% matlab does not support changing a file in the middle, so we need to
% re-write it line by line.

if nargin < 3
   newname = [];
end

[jnk, oldname] = fileparts(oldpath);

fid1 = fopen(oldpath,'r');
fid2 = fopen(newpath,'w');

tline = fgetl(fid1);
done = false;
while ischar(tline)
   tline = strtrim(tline);
   if ~done && ~(isempty(tline) || (length(tline) > 0 && tline(1)=='%') || ...
      (length(tline) > 8 && strcmp(tline(1:8), 'function')))
      write_warning(fid2,oldname, newname);
      done = true;
      %break; % if we break here the rest of the file will not be written
   end
   fprintf(fid2,'%s\n',tline);
   tline = fgetl(fid1);
end
fclose(fid1)
fclose(fid2)

function write_warning(fid, oldname, newname)
if ~isempty(newname)
fprintf(fid,...
   'warning(''EIDORS:deprecated'',''%s is deprecated as of %s. Use %s instead.'');\n\n',...
   date, upper(oldname),upper(newname));
else
   fprintf(fid,...
   'warning(''EIDORS:deprecated'',''%s is deprecated as of %s. '');\n\n',...
   date, upper(oldname));
end

function do_unit_test
fid = fopen('test1.m','w');
fprintf(fid, '\n');
fprintf(fid, ' %% nasty comment \n');
fprintf(fid, 'function test1(x,y)\n');
fprintf(fid, '%% test1(x,y)\n');
fprintf(fid, '\n');
fprintf(fid, ' %% nasty comment \n');
fprintf(fid, '\n');
fprintf(fid, 'code\n');
fprintf(fid, 'code\n');
fprintf(fid, 'code\n');
fclose(fid);

eidors_deprecate('test1','test5');