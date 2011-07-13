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
cd ..
if 0
m2html('mfiles','eidors', 'htmldir','doc','recursive','on',...
    'globalhypertextlinks', 'on','template','blue',...
    'helptocxml', 'on');
else
m2html('mfiles','eidors', 'htmldir','doc','recursive','on',...
    'globalhypertextlinks', 'on','template','frame','index','menu',...
    'helptocxml', 'on');
end
!cp build/doc_template/intro.html doc/
!rsync -r doc htdocs
!rm -rf doc
cd build
!rm -rf m2html
!rm m2html.zip
rmpath([pwd '/m2html']);
