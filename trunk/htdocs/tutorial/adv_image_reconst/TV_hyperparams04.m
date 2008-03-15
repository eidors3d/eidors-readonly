% Generate HTML frame to view $Id: TV_hyperparams04.m,v 1.1 2008-03-15 00:12:20 aadler Exp $

fprintf(fid,'<TABLE><TR><TH>');
for alpha1= alpha1list
   fprintf(fid,'<TH>&alpha;<sub>1</sub>=10<sup>&minus;%3.1f</sup>',alpha1);
end
fprintf(fid,'\n');

for alpha2= alpha2list
   fprintf(fid,'\n<TR>\n<TH>&alpha;<sub>2</sub>=10<sup>&minus;%3.1f</sup>\n',alpha2);
   for alpha1= alpha1list
      fprintf(fid,'<TD><img src="%s">\n', sprintf(name_base,alpha1,alpha2) );
   end
end

fprintf(fid,'\n</TABLE></BODY></HTML>\n');
fclose(fid);


