function elec_idx = elec_rearrange( pattern, newarrange )
% ELEC_REARRANGE: rearrange electrodes for given pattern
%
% Usage: 
%     elec_idx = elec_rearrange( pattern, newarrange );
%     fmdl.electrode = fmdl.electrode( elec_idx );
%
% pattern    = [Elecs/row, rows]  (eg [16,2])
%   or a structure with (p.n_elecs, p.n_layers)
%
% newarrange = 'square' 
%              'zigzag' / 'odd_even' 
%              'none' 
%
%
% (C) 2017 Andy Adler. License: GPL v2 or v3. $Id$

   if ischar(pattern) && strcmp(pattern,'UNIT_TEST'); do_unit_test; return; end

   if isnumeric(pattern)
      p.n_elecs = pattern(1);
      p.n_layers= pattern(2);
   elseif isstruct(pattern)
      p = pattern;
   else
      error('Can''t interpret "pattern" input');
   end

   if p.n_layers >= 3;
      error('elec_rearrange: can''t interpret n_layers >=2');
   end

   idx = reshape(1:(p.n_elecs*p.n_layers),p.n_layers,p.n_elecs);
   switch newarrange
     case 'none';      return
     case 'odd_even';
        idx = reshape(idx,2,[])';
     case 'square'  ;
        idx = reshape(idx,2,[])';
        idx(2:2:end,:) = fliplr(idx(2:2:end,:));
     otherwise;
          error('case %s not understood',newarrange);
   end
  elec_idx = idx(:);

function do_unit_test
   sqV = [1,2;4,3;5,6;8,7;9,10;12,11;13,14;16,15]';


   elec= reshape([1:16],8,2)';
   idx = elec_rearrange(size(elec'),'square');
   idx = reshape(idx,size(elec'))';
   unit_test_cmp('Sq [8,2]', idx, sqV);

   p.n_elecs = 8; p.n_layers = 2;
   idx = elec_rearrange(p,'square');
   idx = reshape(idx,size(elec'))';
   unit_test_cmp('Sq [8,2]', idx, sqV);

   elec= [1:32]';
   idx = elec_rearrange(size(elec),'square');
   idx= reshape(idx,[16,2])';
   unit_test_cmp('Sq [16,2]', idx(:,1:8), sqV);

   imdl = mk_common_model('n3r2',[16,2]); fmdl = imdl.fwd_model;
   idx = elec_rearrange([16,2],'square');
   fmdl.electrode(idx) = fmdl.electrode;
   show_fem_enhanced(fmdl,[0,1]);
   view(10,82);
   
