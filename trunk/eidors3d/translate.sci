files= unix_g("ls *.m");
nf = size(files,1);
for k = 1:nf, mfile2sci(files(k),"scifiles",%f,%f),end