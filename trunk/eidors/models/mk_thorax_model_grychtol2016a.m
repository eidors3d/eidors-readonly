function img = mk_thorax_model_grychtol2016a(stimpat)
% Builds the male thorax models from Grychtol, Müller and Adler (2016)
% Accepts a string parameter:
%    '2x16_planar'   - 2 rings of 16 electrodes, planar stimulation
%    '2x16_odd-even' - 2 rings, each stimulation is cross-plane 
%    '2x16_square'   - 2 rings, every second stimulation is cross-plane
%    '2x16_adjacent' - like square, but with adjacent stimulation
%    '1x32_ring'     - single ring of 32 electrodes
%
% CITATION_REQUEST:
% AUTHOR: Bartlomiej Grychtol, Beat Müller and Andy Adler
% TITLE: 3D EIT image reconstruction with GREIT
% JOURNAL: Physiological Measurement
% VOL: 37
% YEAR: 2016

% (C) 2016 Bartlomiej Grychtol. License: GPL version 2 or 3
% $Id$

citeme(mfilename);

if nargin<1 || isempty(stimpat)
   stimpat = '1x32_ring';
end
   

eth16 = 360*cumsum([0.2 0.4*ones(1,7) 0.5 0.4*ones(1,7)])/6.5 - 90; eth16 = eth16';
eth32 = 360*cumsum([0.1 0.2*ones(1,15) 0.25 0.2*ones(1,15)])/6.45 - 90; eth32 = eth32';
ep = eth16; ep(:,2) = 150;
ep(17:48,1) = eth32; ep(17:48,2) = 175;
ep(49:64,1) = eth16; ep(49:64,2) = 200;
mdl = mk_thorax_model('male',ep,[5 0 .5],10);
mdl.name = sprintf(['Thorax mesh from Grychtol, Müller and Adler (2016) '... 
                     '- %s electrode config'],stimpat);

emap = get_elec_map(stimpat);
mdl.electrode = mdl.electrode(emap);
mdl.stimulation = get_stim_pattern(stimpat);


opt = organ_options;
img = mk_lung_image(mdl,opt);
img.fwd_model.normalize_measurements = false;

img.name = sprintf(['Thorax model from Grychtol, Müller and Adler (2016) '... 
                     '- %s electrode config'],stimpat);

end


function  opt = organ_options
   opt.bkgnd_val = 1  ;
   opt.lung_val =  .2;
   opt.heart_val = 1.5;
   opt.left_lung_center =  [ 75 25 100];
   opt.right_lung_center = [-75 25 100];
   opt.left_lung_axes = [80 100 250];
   opt.right_lung_axes = [80 100 250];
   opt.heart_center = [20 -25 200];
   opt.heart_axes = [50 60 75];
   opt.diaphragm_center = [0 -50 0];
   opt.diaphragm_axes = [220 190 120];
end

function stim = get_stim_pattern(str)
   switch(str)
      case '2x16_planar'
         stim = mk_stim_patterns(32,1,[0 6],[0 6],{'no_meas_current','no_rotate_meas'},1);
      case {'2x16_odd-even', '2x16_square', '1x32_ring'}
         stim = mk_stim_patterns(32,1,[0 5],[0 5],{'no_meas_current','no_rotate_meas'},1);
      case '2x16_adjacent'
         stim = mk_stim_patterns(32,1,[0 1],[0 1],{'no_meas_current','no_rotate_meas'},1);  
      otherwise
         error('Stim pattern string not understood. Available strings are: \n%s', ...
            sprintf('%s\n', pattern_list));
   end
end

function ls = pattern_list
   ls = {
      '2x16_planar'
      '2x16_odd-even'
      '2x16_square'
      '2x16_adjacent'
      '1x32_ring'
      };
end

function map = get_elec_map(str)

   switch str
      case {'2x16_odd-even', '2x16_planar'}
         map = oddeven32;
      case {'2x16_square', '2x16_adjacent'}
         map = square32;
      case '1x32_ring'
         map = ring32;
      otherwise
         error('No such electrode map');
   end
end

function m = square32
   o = [48 1 -48 1];
   o = repmat(o,1,8);
   m = zeros(1,32);
   m(1) = 1;
   for i = 2:32
      m(i) = m(i-1) + o(i-1);
   end
end


function m = oddeven32
   odd = 49:64; % top layer
   even = 1:16; % bottom layer
   m = [odd; even];
   m = m(:)';
end

function m = ring32
   m = 17:48;
end