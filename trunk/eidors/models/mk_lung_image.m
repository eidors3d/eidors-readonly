function img = mk_lung_image(mdl, opt, cmd)
%MK_LUNG_IMAGE create an image with lung and heart contrasts from a model
% MK_LUNG_IMAGE defines the lungs, heart and diaphragm as ellipsoids and
% creates an image with these organs out of the supplied eidors fwd_model.
%
% Typical usage:
%
% IMG = MK_LUNG_IMAGE(MDL) uses default options for shape, position and
%  contrast. 
% IMG = MK_LUNG_IMAGE(MDL, OPT) uses the specified options (see below).
% IMG = MK_LUNG_IMAGE(MDL, OPT, 'show') also generate a figure with 4 views
%  of the produced image
%
% Additional utilities:
%
% OPT = MK_LUNG_IMAGE(MDL,'options') returns the default options for the 
%  given model, and does not generate an image.
% MK_LUNG_IMAGE(MDL, OPT, 'diagnostic') shows the supplied model together 
%  with the organ ellipsoids as defined by the supplied options. This usage
%  has no return parameters.
%
% The options input OPT is a struct with the following fields and defaults:
%             bkgnd_val: 1          - image value for the background
%              lung_val: 0.2000     - image value for the lungs
%             heart_val: 1.5000     - image value for the hearts
%      left_lung_center: ctr + [  .24    0   -1/3]*scale
%     right_lung_center: ctr + [ -.24    0   -1/3]*scale
%          heart_center: ctr + [ 0.07  -1/6     0]*scale
%      diaphragm_center: ctr + [ 0     -1/4  -2/3]*scale
%            heart_axes: [ 1/6  1/5   1/4 ]*scale
%        left_lung_axes: [ 8/30 1/3  25/30]*scale
%       right_lung_axes: [ 8/30 1/3  25/30]*scale
%        diaphragm_axes: [22/30 19/30  2/5]*scale
% where ctr is the center of the bounding box of the model, and scale is 
% the smaller of the bounding box's x dimension divided by 7/6 and its y 
% dimension. Each organ is an ellipsoid with the corresponding center and
% axes. 
% NOTE that only options that deviate from defaults need be specified.
%
% Example:
%    mdl = mk_thorax_model('female');       % female thorax
%    img = mk_lung_image(mdl);
%    subplot(131)
%    show_fem(img) 
%    % it seems that the heart definition goes outside the model
%    % let's correct the default options
%    opt = mk_lung_image(mdl, 'options');   % get default options 
%    corr = 20;
%    opt.heart_center(2)        = opt.heart_center(2)      + corr;
%    opt.left_lung_center(2)    = opt.left_lung_center(2)  + corr;
%    opt.right_lung_center(2)   = opt.right_lung_center(2) + corr;
%    opt.diaphragm_center(2)    = opt.diaphragm_center(2)  + corr;
%    img = mk_lung_image(mdl, opt);
%    subplot(132)
%    show_fem(img)
%    % it's hard to tell because the model has very big elements.
%    % let's make a diagnostic plot to verify the organ definitions fit inside 
%    % the model (independent of model maxh)
%    subplot(133)
%    mk_lung_image(mdl,opt,'diagnostic');  
%
% CITATION_REQUEST:
% AUTHOR: Bartlomiej Grychtol, Beat Müller and Andy Adler
% TITLE: 3D EIT image reconstruction with GREIT
% JOURNAL: Physiological Measurement
% VOL: 37
% YEAR: 2016
%
% See also: ELEM_SELECT, MK_THORAX_MODEL

% (C) 2015-2016 Bartlomiej Grychtol. License: GPL version 2 or 3
% $Id$

citeme(mfilename);

if nargin>0 && ischar(mdl) && strcmp(mdl,'UNIT_TEST')
   do_unit_test; return;
end

if nargin == 1, opt = struct; end
returnopt = false;
if nargin == 2 && ischar(opt) && strcmp(opt,'options')
   returnopt = true;
   clear opt;
   opt = struct;
end
opt = parse_opt(mdl,opt);
if returnopt
   img = opt;
   return
end

if nargin == 3 && ischar(cmd) && strcmp(cmd, 'diagnostic')
   show_organs(mdl,opt);
   return
end

elem_data = define_elem_data(mdl, opt);

img = mk_image(mdl,elem_data);

if nargin==3 && ischar(cmd) && strcmp(cmd, 'show')
   show4views(img);
end
end

function elem_data = define_elem_data(mdl, opt)

   el_str = '(x-%f).^2/%f^2 + (y-%f).^2/%f^2 + (z-%f).^2/%f^2 <= 1';
   pack = @(x,y) mat2cell(reshape([x' y']',1,[])',ones(1, 2*numel(x)))';


   args = pack(opt.left_lung_center, opt.left_lung_axes);
   ll_str = sprintf(el_str, args{:});
   ll_fcn = inline(ll_str,'x','y','z');
   left_lung = elem_select(mdl,ll_fcn);

   args = pack(opt.right_lung_center, opt.right_lung_axes);
   rl_str = sprintf(el_str, args{:});
   rl_fcn = inline(rl_str,'x','y','z');
   right_lung = elem_select(mdl,rl_fcn);

   lung = min(1,left_lung + right_lung);

   args = pack(opt.diaphragm_center, opt.diaphragm_axes);
   dia_str = sprintf(el_str, args{:});
   dia_fcn = inline(dia_str,'x','y','z');
   diaphragm = elem_select(mdl,dia_fcn);

   lung = max(0,lung - diaphragm);

   args = pack(opt.heart_center, opt.heart_axes);
   hrt_str = sprintf(el_str, args{:});
   hrt_fcn = inline(hrt_str,'x','y','z');
   heart = elem_select(mdl, hrt_fcn);

   lung = max(0, lung - heart);

   elem_data = ...
       opt.bkgnd_val ...
       - (opt.bkgnd_val-opt.lung_val)*lung ...
       - (opt.bkgnd_val-opt.heart_val)*heart;
end

function show_organs(mdl, opt)
   show_fem(mdl)
   hold on
   pack = @(x) mat2cell(x',ones(1,numel(x)))';
   args = pack([opt.left_lung_center opt.left_lung_axes]);
   [x y z] = ellipsoid(args{:});
   surf(x,y,z,'FaceColor','b');

   args = pack([opt.right_lung_center opt.right_lung_axes]);
   [x y z] = ellipsoid(args{:});
   surf(x,y,z,'FaceColor','b');

   args = pack([opt.heart_center opt.heart_axes]);
   [x y z] = ellipsoid(args{:});
   surf(x,y,z,'FaceColor','r');


   args = pack([opt.diaphragm_center opt.diaphragm_axes]);
   [x y z] = ellipsoid(args{:});
   surf(x,y,z,'FaceColor','g');
   hold off
end

function show4views(img)
    subplot(221)
    h = show_fem(img);
    set(h,'linewidth',.1)
    subplot(222)
    h = show_fem(img);
    set(h,'linewidth',.1)
    view([0 90])
    subplot(223)
    h = show_fem(img);
    set(h,'linewidth',.1)
    view([90,0])
    subplot(224)
    h = show_fem(img);
    set(h,'linewidth',.1)
    view([0 0])
end

function opt = parse_opt(fmdl,opt)
   minnode = min(fmdl.nodes);
   maxnode = max(fmdl.nodes);
   bb = maxnode-minnode;
   scale = min( bb(1:2) ./ [7/6 1]); % the default values below fit a model of these proportions best
   ctr = minnode + bb/2;

   opt = set_default(opt, 'bkgnd_val',          1  );
   opt = set_default(opt, 'lung_val' ,           .2);
   opt = set_default(opt, 'heart_val',          1.5);
    
   opt = set_default(opt, 'left_lung_center',  ctr + [  .24    0   -1/3]*scale);
   opt = set_default(opt, 'right_lung_center', ctr + [ -.24    0   -1/3]*scale);
   opt = set_default(opt, 'heart_center',      ctr + [ 0.07  -1/6     0]*scale);
   opt = set_default(opt, 'diaphragm_center',  ctr + [ 0     -1/4  -2/3]*scale);
   
   opt = set_default(opt, 'heart_axes',        [ 1/6  1/5   1/4 ]*scale);
   opt = set_default(opt, 'left_lung_axes' ,   [ 8/30 1/3  25/30]*scale);
   opt = set_default(opt, 'right_lung_axes',   [ 8/30 1/3  25/30]*scale);
   opt = set_default(opt, 'diaphragm_axes',    [22/30 19/30  2/5]*scale);
end

function s = set_default(s, fname, defval)
   if ~isfield(s,fname)
       s.(fname) = defval;
   end
end


function do_unit_test
mdl = ng_mk_ellip_models([400 180 150 10],4,10);
img = mk_lung_image(mdl);
subplot(231)
show_fem(img);

% sclaing should have no effect
mdl.nodes = 3*mdl.nodes;
img = mk_lung_image(mdl);
subplot(232)
show_fem(img);

% shorter model should show a "slice"
mdl = ng_mk_ellip_models([200 180 150 10],4,10);
img = mk_lung_image(mdl);
subplot(233)
show_fem(img);

% it should still fit into the male torso
mdl = mk_thorax_model('male');
img = mk_lung_image(mdl);
subplot(234)
show_fem(img);

% the female model seems to need adjustements
mdl = mk_thorax_model('female');       % female thorax
opt = mk_lung_image(mdl, 'options');   % 
corr = 20;
opt.heart_center(2)        = opt.heart_center(2)      + corr;
opt.left_lung_center(2)    = opt.left_lung_center(2)  + corr;
opt.right_lung_center(2)   = opt.right_lung_center(2) + corr;
opt.diaphragm_center(2)    = opt.diaphragm_center(2)  + corr;
subplot(235)
% make a diagnostic plot to verify the organ definitions fit inside the
% model (independent of model maxh)
mk_lung_image(mdl,opt,'diagnostic');  
% show_fem(img);

return

mdl = mk_thorax_model('grychtol2016a_male_thorax');
img = mk_lung_image(mdl);
subplot(236)
show_fem(img);

end