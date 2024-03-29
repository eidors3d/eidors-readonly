% Generate HTML frame to view $Id$

fid= fopen('TV-params-NSR=0.html','w');

a=sprintf('%calpha;',38); % alpha
m=sprintf('%cminus;',38); % alpha
s=sprintf('%c',60); % less than 
e=sprintf('%c',62); % greater than
tr= [s,'TR',e]; etr= [s,'/TR',e];
th= [s,'TH',e]; eth= [s,'/TH',e];
td= [s,'TD',e]; etd= [s,'/TD',e];
sub=[s,'SUB',e];esub=[s,'/SUB',e];
sup=[s,'SUP',e];esup=[s,'/SUP',e];


fprintf(fid,[s,'TABLE',e,tr,th]);
for alpha1= alpha1list
   fprintf(fid,[th,a,sub,'1',esub,'=10',sup,m,'%3.1f',esup],alpha1);
end
fprintf(fid,'\n');

for alpha2= alpha2list
   fprintf(fid,['\n',tr,'\n',th,a,sub,'2',esub,'=10', ...
                sup,m,'%3.1f',esup,'\n'],alpha2);
   for alpha1= alpha1list
      fprintf(fid,[td,s,'img src="%s"',e,'\n'], sprintf(name_base,alpha1,alpha2) );
   end
end

fprintf(fid,['\n',s,'/TABLE',e,'\n']);
fclose(fid);
