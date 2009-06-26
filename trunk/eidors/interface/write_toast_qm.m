function write_toast_qm(fname, fmdl)
fid = fopen(fname,'w');
dims = size(fmdl.nodes,2);
elecs= length(fmdl.electrode);
elecpos= zeros(elecs,dims);
for i=1:elecs
    elnodes = fmdl.electrode(i).nodes;
    elecpos(i,:) = mean( fmdl.nodes(elnodes,:),1);
end

fprintf(fid,'QM file %dD\n',dims);
fprintf(fid,'Dimension %d\n',dims);
fprintf(fid,'SourceList %d fixed\n',elecs);
fprintf(fid,'%6.4g     %6.4g\n', elecpos');
fprintf(fid,'MeasurementList %d\n',elecs);
fprintf(fid,'%6.4g     %6.4g\n', elecpos');
fprintf(fid,'LinkList\n');
for i=1:elecs
  fprintf(fid,'%d:',elecs);
  fprintf(fid,'%d ',0:elecs-1);
  fprintf(fid,'\n');
end

fclose(fid);