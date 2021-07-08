function opt = ng_read_opt(fname)
%NG_READ_OPT Read Netgen's ng.opt 
%  NG_READ_OPT, without inputs, reads ng.opt in current directory
%
%  NG_WRITE_OPT(PATH) reads the file specified in PATH. IF PATH is a 
%  directory, looks for ng.opt.
% 
%  See also NG_WRITE_OPT

% (C) 2021 Bartlomiej Grychtol. License: GPL version 2 or version 3
% $Id$

% if input is 'UNIT_TEST', run tests
if nargin == 1 && ischar(fname) && strcmp(fname,'UNIT_TEST') 
   do_unit_test; return; end

if nargin == 0 
    fname = 'ng.opt';
end

if isfolder(fname)
    fname = [fname filesep 'ng.opt'];
end

opt = struct();

fid = fopen(fname,'r'); 
tline = fgetl(fid);
while ischar(tline)
    [key, val] = strtok(tline,' ');
    strval = strip(val);  
    numval = str2double(strval);
    if isnan(numval)
        eval(sprintf('opt.%s = ''%s'';',key,strval));
    else
        eval(sprintf('opt.%s = %f;',key,numval));
    end
    tline = fgetl(fid);
end
fclose(fid);



function do_unit_test
    opt = ng_write_opt();
    opt.meshoptions.fineness = 4;
    ng_write_opt(opt);
    opt = ng_read_opt();
    unit_test_cmp('simple test',opt.meshoptions.fineness, 4);
    
    path = cd();
    opt = ng_read_opt(path);
    unit_test_cmp('dir test',opt.meshoptions.fineness, 4);
    
    path = cd();
    opt = ng_read_opt([path filesep 'ng.opt']);
    unit_test_cmp('path test',opt.meshoptions.fineness, 4);
    
    
    
