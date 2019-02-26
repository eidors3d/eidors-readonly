function emphasize_electrode_text(elecs);
% EMPHASIZE_ELECTRODE_TEXT: ephazise the electrode
%   text on a "show_fem" plot
% emphasize_electrode_text(elecs)
%   elecs => which elecs to emphasize (default, all)
%
% (C) 2005-2018 Andy Adler. License: GPL version 2 or version 3
% $Id$

hh = get(gca,'Children');

bgcolour = [1,1,1];
tmp = inv(view);
dirn = tmp(1:3,3); dirn = -dirn'/norm(dirn);

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
  set(hi, 'Position', pn + dirn*10);
end

