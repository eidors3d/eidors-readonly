function mdl = ng_mk_2d_models(shape)


if ischar(shape) && strcmp(shape, 'UNIT_TEST'), do_unit_test, return, end 
if ~iscell(shape)
   shape = {shape}
end

points = [];
for i = 1:length(shape)
   lp = length(points);
   ls = length(shape{i});
   points = [points; shape{i}];
   seg{i} = repmat([0 1],ls,1) + lp + repmat((1:ls)',1,2);
   seg{i}(end,2) = lp + 1;
end
write_in2d_file('tmp.in2d',points, seg);


function write_in2d_file(fname,points, seg)
fid = fopen(fname,'w');
fprintf(fid, '%s\n','splinecurves2dv2');
fprintf(fid, '%d\n',6); % global grading factor, 6 should force use of ng.opt
fprintf(fid, '%s\n','points');
for i = 1:length(points)
   fprintf(fid, '%d   %f   %f\n',i,points(i,:));
end
fprintf(fid,'%s\n','segments');
% here we assume the first loop is the boundary, all the others are holes
domains = [ 1 0];
for i = 1:length(seg)
   if i > 1
      domains = [0 1];
   end
   for j = 1:length(seg{i})
      fprintf(fid,'%d   %d   %d   %d   %d -bc=%d\n',domains, 2, seg{i}(j,:),i);
   end
end
fclose(fid);

function do_unit_test
xy = [0 0;  1 0; 1 1; 0 1];
ng_mk_2d_models({xy, 0.25 + 0.5*xy});