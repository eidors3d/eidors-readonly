function a = common_colourbar(ha,varargin)
%COMMON_COLOURBAR make a joint colourbar for several subplots
% a = common_colourbar(ha,img) adds a common colourbar to the subplots
% identified in the array of handles ha, passing the img input to
% eidors_colourbar for actual rendering, and returning a handle to a new
% parent axis of the colourbar (a).
%
% a = common_colourbar(ha, ...) uses all inputs beside ha to call
% eidors_colourbar.
%
% common_colourbar('clear',a) removes the colourbar previously added by
% a call to common_colourbar.

% (C) 2013 Bartlomiej Grychtol. License: GPL version 2 or 3.
% $Id$

if ischar(ha)
   switch ha
      case 'UNIT_TEST'
         do_unit_test;
         return;
      case 'clear'
         img = varargin{1};
         delete(img);
      otherwise
         error('huh?');
   end
end
   
fig = get(ha(1),'Parent');
ca = get(fig,'CurrentAxes');
%1. Figure out position limits
x_min = inf; x_max = -inf;
y_min = inf; y_max = -inf;
for i = 1:length(ha)
   u = get(ha(i),'Units');
   set(ha(i),'Units','points')
   pos = get(ha(i), 'Position');
   x_min = min(x_min, pos(1));
   x_max = max(x_max, pos(1)+pos(3));
   y_min = min(y_min, pos(2));
   y_max = max(y_max, pos(2)+pos(4));
   set(ha(i),'Units',u);
end
%2. Create an axis that overlaps them all
pos = [x_min y_min x_max-x_min y_max-y_min];
% this is how Matlab R2011a calculates space for the colorbar
pp = get(gca,'Position');
sz = min(max(pp(3:4)*0.3,20),40);
% so we pre-emptively extend the axis by this amount
pos(3) = pos(3) + sz(1);

a = axes('Units','points','Position',pos);
axis off

%3. Add eidors colorbar
eidors_colourbar(varargin{:});

%4. Minimize impact on the figure from user's perspective

% put the new axis behind everything else
set(fig,'Children',circshift(get(fig,'Children'),-1));

% return focus to the previous axes
set(fig,'CurrentAxes',ca)




function do_unit_test
clf
for i = 1:6
h(i)  = subplot(2,3,i);
end
common_colourbar(h([1 2 4 5]),2,0);
common_colourbar(h(4:6),5,1);
eidors_msg('CHECK GRAPHICS',0);
