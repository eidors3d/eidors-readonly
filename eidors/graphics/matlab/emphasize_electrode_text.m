function emphasize_electrode_text;
% EMPHASIZE_ELECTRODE_TEXT: ephazise the electrode
%   text on a "show_fem" plot
% emphasize_electrode_text(hh)
%
% (C) 2005-2018 Andy Adler. License: GPL version 2 or version 3
% $Id$

hh = get(gca,'Children');

bgcolour = [1,1,1];
tmp = view;
dirn = tmp(1:3,3); dirn = -dirn'/norm(dirn) 

for hi = hh(:)'
  if ~strcmp(get(hi,'Type'),'text'); 
     continue;
  end
  set(hi, 'BackgroundColor', bgcolour);
  set(hi, 'FontSize', 12);
  pn = get(hi, 'Position');
  set(hi, 'Position', pn + dirn);
end

