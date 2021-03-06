% run in the build directory
%STUPID )(*&()*&)#@$ can't work with https
%urlwrite('https://www.artefact.tk/software/matlab/m2html/m2html.zip',...
%    'm2html.zip');
!wget --no-check-certificate https://www.artefact.tk/software/matlab/m2html/m2html.zip
unzip('m2html.zip');
addpath([pwd '/m2html']);
tpl = 'blue';
try 
!chmod u+w -R m2html
end
copyfile('doc_template/index.html','m2html/templates/frame/index.html','f')
copyfile('doc_template/matlabicon.gif','m2html/templates/frame/matlabicon.gif','f')
copyfile('doc_template/matlabicon.gif','m2html/templates/blue/matlabicon.gif','f')
copyfile('doc_template/mfile.tpl','m2html/templates/frame/mfile.tpl','f')

VERSION = 0; % 1 for MATLAB docs, 0 for SOURCEFORGE


 cd ..
% FOR THE BUILD PROCESS, AFTER RUNNING MAKE STEPS 0 .. 9
if VERSION
m2html('mfiles','eidors', 'htmldir','doc','recursive','on',...
    'globalhypertextlinks', 'on','template','blue',...
    'helptocxml', 'on');
else
copyfile('build/doc_template/m2html.m','build/m2html','f');
m2html('mfiles','eidors', 'htmldir','doc','recursive','on',...
    'globalhypertextlinks', 'on','template','frame','index','menu',...
    'helptocxml', 'on');
 copyfile('build/doc_template/intro.html','doc/intro.html','f')
end
!rsync -r doc htdocs
!rm -rf doc
if VERSION
   cd htdocs/doc
   p = cd;
   %eidors must be on the path before this can run
   builddocsearchdb(p);
   cd ../..
end
cd build
!rm -rf m2html
!rm m2html.zip
rmpath([pwd '/m2html']);

% If building for the release,  move the created doc directory to the release
