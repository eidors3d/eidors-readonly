urlwrite('http://www.artefact.tk/software/matlab/m2html/m2html.zip',...
    'm2html.zip');
!unzip m2html
addpath([pwd '/m2html']);
tpl = 'blue';
!chmod u+w -R m2html
!cp doc_template/index.html m2html/templates/frame
!cp doc_template/matlabicon.gif m2html/templates/frame
!cp doc_template/matlabicon.gif m2html/templates/blue
if 1
m2html('mfiles','eidors', 'htmldir','doc','recursive','on',...
    'globalhypertextlinks', 'on','template','blue',...
    'helptocxml', 'on');
else
m2html('mfiles','eidors', 'htmldir','doc','recursive','on',...
    'globalhypertextlinks', 'on','template','frame','index','menu',...
    'helptocxml', 'on');
end
!rsync -r doc htdocs
!rm -rf doc
!rm -rf m2html
!rm m2html.zip
rmpath([pwd '/m2html']);