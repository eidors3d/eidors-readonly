function emphasize_electrode_text(elecs);
% EMPHASIZE_ELECTRODE_TEXT: ephazise the electrode
%   text on a "show_fem" plot
% emphasize_electrode_text(elecs)
%   elecs => which elecs to emphasize (default, all)
%
% (C) 2005-2018 Andy Adler. License: GPL version 2 or version 3
% $Id$

% FIXME: the view matrix is not easy to interpret. Fix dirn

if (nargin>=1) && ischar(elecs) && strcmp(elecs,'UNIT_TEST'); do_unit_test; return; end

hh = get(gca,'Children'); axis vis3d

bgcolour = [1,1,1];
% Don't really understand the view matrix
%tmp = inv(view);
 tmp = view;
dirn = [0,0,1,0]*tmp; % move into view's z dir
dirn = dirn(1:3);
dirn = -dirn(:)'/norm(dirn);

if nargin==0;
  elecs= 1:length(hh)'; % all of them
end

k=0;for hi = hh(:)'; k=k+1;
  if ~strcmp(get(hi,'Type'),'text'); 
     continue;
  end
  if ~any(str2num(get(hi,'String'))==elecs);
      continue;
  end
  set(hi, 'BackgroundColor', bgcolour);
  set(hi, 'FontSize', 12);
  pn = get(hi, 'Position');
  set(hi, 'Position', pn + dirn*1 );
end

function do_unit_test
  fmdl = getfield(mk_common_model('n3r2',[16,2]),'fwd_model');
  subplot(221);
  show_fem(fmdl,[0,1]); view(0,0); drawnow

  emphasize_electrode_text();
  subplot(222);
  show_fem(fmdl,[0,1]); view(0,80); drawnow

  emphasize_electrode_text();
  subplot(223);
  show_fem(fmdl,[0,1]); view(50,80); drawnow

  emphasize_electrode_text();
