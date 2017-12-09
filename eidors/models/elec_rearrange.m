function [elec_idx, new_fmdl] = elec_rearrange( pattern, newarrange, fwd_model )
% ELEC_REARRANGE: rearrange electrodes for given pattern
%
% Usage: 
%     elec_idx = elec_rearrange( pattern, newarrange );
%     fmdl.electrode( elec_idx ) = fmdl.electrode;
%   OR;
%     [~,new_fmdl] = elec_rearrange( pattern, newarrange, fwd_model );
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
   switch lower(newarrange)
     case 'none'; % nothing required
     case {'odd_even','zigzag'};
        idx = reshape(idx,2,[])';
     case 'square'  ;
        idx = reshape(idx,2,[])';
        idx(2:2:end,:) = fliplr(idx(2:2:end,:));
     otherwise;
          error('case %s not understood',newarrange);
   end
  elec_idx = idx(:);

  if nargin >= 3
     new_fmdl = fwd_model;
     new_fmdl.electrode(elec_idx) = new_fmdl.electrode;
  end

function do_unit_test
   sqV = [1,2;4,3;5,6;8,7;9,10;12,11;13,14;16,15]';
   oeV = [1,2;3,4;5,6;7,8;9,10;11,12;13,14;15,16]';


   elec= reshape([1:16],8,2)';
   idx = elec_rearrange(size(elec'),'square');
   unit_test_cmp('Sq#1 [8,2]', reshape(idx,size(elec'))', sqV);

   idx = elec_rearrange(size(elec'),'odd_even');
   unit_test_cmp('OE#1 [8,2]', reshape(idx,size(elec'))', oeV);

   idx = elec_rearrange(size(elec'),'zigZaG');
   unit_test_cmp('OE#2 [8,2]', reshape(idx,size(elec'))', oeV);

   idx = elec_rearrange(size(elec'),'none');
   unit_test_cmp('OE#2 [8,2]', idx', 1:16);

   p.n_elecs = 8; p.n_layers = 2;
   idx = elec_rearrange(p,'square');
   idx = reshape(idx,size(elec'))';
   unit_test_cmp('Sq#2 [8,2]', idx, sqV);

   elec= [1:32]';
   idx = elec_rearrange(size(elec),'square');
   idx= reshape(idx,[16,2])';
   unit_test_cmp('Sq [16,2]', idx(:,1:8), sqV);

   imdl = mk_common_model('n3r2',[16,2]); fmdl1= imdl.fwd_model;
   subplot(221);
   idx = elec_rearrange([16,2],'square');
   fmdl1.electrode(idx) = fmdl1.electrode;
   show_fem_enhanced(fmdl1,[0,1]);
   view(10,82);

   subplot(222);
   [~,fmdl2] = elec_rearrange([16,2],'square', imdl.fwd_model);
   show_fem_enhanced(fmdl2,[0,1]);
   view(10,82);

   unit_test_cmp('fmdl#1', [fmdl1.electrode([5,9]).nodes], [156 157 219 220 128 129 191 192]);
   unit_test_cmp('fmdl#2', fmdl1.electrode, fmdl2.electrode) 
