function eidors_rename(oldname, newname)
% eidors_rename(oldname, newname)
% Renames the function oldname to newname and replaces all calls in eidors.
% Ignores files in the depricated directory.

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

curdir = cd;
cd(olddir);
system(sprintf('%s svn mv %s.m %s.m', LDP, oldname, newname));
cd(eidors_dir);
match = '-path ./deprecated -prune -o -path ./tools -prune -o -type f -iname "*.m"'; 
system(sprintf(...
   'find . %s -print0 | xargs -0 sed  -i "s/%s/%s/g;s/%s/%s/g;"',...
   match, oldname, newname,upper(oldname),upper(newname)));
match =  '-type f -iname "*.m"';
system(sprintf(...
   'find ../dev %s -print0 | xargs -0 sed  -i "s/%s/%s/g;s/%s/%s/g;"',...
   match, oldname, newname,upper(oldname),upper(newname)));
system(sprintf(...
   'find ../htdocs/tutorial %s -print0 | xargs -0 sed  -i "s/%s/%s/g;s/%s/%s/g;"',...
   match, oldname, newname,upper(oldname),upper(newname)));
cd(curdir)

function ldp = LDP;
   ldp = 'LD_LIBRARY_PATH="" ';

function do_unit_test
fid = fopen('test5.m','w');
fprintf(fid, 'function test5(x,y)\n');
fprintf(fid, '%% test5(x,y)\n');
fclose(fid);

fid = fopen('test2.m','w');
fprintf(fid, 'function test2(x,y)\n');
fprintf(fid, '%%test function 2\n');
fprintf(fid, 'test2(x,y);');
fclose(fid);

eidors_rename('test5', 'test3');
