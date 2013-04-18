% run in the build directory
urlwrite('http://www.artefact.tk/software/matlab/m2html/m2html.zip',...
    'm2html.zip');
unzip('m2html.zip');
addpath([pwd '/m2html']);
tpl = 'blue';
!chmod u+w -R m2html
!cp doc_template/index.html m2html/templates/frame
!cp doc_template/matlabicon.gif m2html/templates/frame
!cp doc_template/matlabicon.gif m2html/templates/blue
!cp doc_template/mfile.tpl m2html/templates/frame

VERSION = 1; % 1 for MATLAB docs, 0 for SOURCEFORGE


 cd ..
% FOR THE BUILD PROCESS, AFTER RUNNING MAKE STEPS 0 .. 9
%cd ~/eidors-release/eidors
if VERSION
m2html('mfiles','eidors', 'htmldir','doc','recursive','on',...
    'globalhypertextlinks', 'on','template','blue',...
    'helptocxml', 'on');
else
!cp build/doc_template/m2html.m m2html
m2html('mfiles','eidors', 'htmldir','doc','recursive','on',...
    'globalhypertextlinks', 'on','template','frame','index','menu',...
    'helptocxml', 'on');
 !cp build/doc_template/intro.html doc/
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
