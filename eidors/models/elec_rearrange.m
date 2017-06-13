function elec_idx = elec_rearrange( pattern, newarrange )
% ELEC_REARRANGE: rearrange electrodes for given pattern
%
% Usage: 
%     elec_idx = elec_rearrange( pattern, newarrange );
%     fmdl.electrode = fmdl.electrode( elec_idx );
%
% pattern    = [Elecs/row, rows]  (eg [16,2])
%   or a structure with (p.n_elecs, p.layers)
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
      p.layers  = pattern(2);
   elseif isstruct(pattern)
      p = pattern;
   else
      error('Can''t interpret "pattern" input');
   end

   if p.layers >= 3;
      error('elec_rearrange: can''t interpret p.layers >=2');
   end

   idx = 1:(p.n_elecs*p.layers);
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
   idx = elec_rearrange([8,2],'square');
   elec= [1:16]; elec = elec(idx); elec= reshape(elec,[8,2]);
   sqV = [1,2;4,3;5,6;8,7;9,10;12,11;13,14;16,15];
   unit_test_cmp('Sq [8,2]', elec, sqV);

   idx = elec_rearrange([16,2],'square');
   elec= [1:32]; elec = elec(idx); elec= reshape(elec,[16,2]);
   unit_test_cmp('Sq [16,2]', elec(1:8,:), sqV);
