function cb_axis = add_colorbar(ax,c1,c2,title,varargin)
%ax - parent axis

% usage add_colorbar(gca, [0 1 0],[1 0 1],'title','LineWidth',3)

gap = 1;
cb_width_pt = 10; %points


h_fig = get(ax,'Parent');
org_unit = get(h_fig,'Units');
set(h_fig,'Units','Points');
fig_pos = get(h_fig,'Position');
set(h_fig,'Units',org_unit);

cb_width = cb_width_pt/fig_pos(3);
ax_pos = get(ax,'Position');
old_ax_width = ax_pos(3);
ax_pos(3) = ax_pos(3) - (1+gap) * cb_width;
set(ax,'Position',ax_pos);


ud = get(ax,'UserData');
if ~isempty(ud)
   for i = 1:length(ud);
       pos = get(ud(i),'Position');
       pos(4) = ax_pos(4);
       set(ud(i),'Position',pos);
   end
end


cb_pos(1) = ax_pos(1) + ax_pos(3) + gap*cb_width;
cb_pos(2) = ax_pos(2);
cb_pos(3) = cb_width;
% if strcmp(get(ax, 'DataAspectRatio'), 'auto')
     cb_pos(4) = ax_pos(4);
% else
%     cb_pos(4) = ax_pos(4)* ax_pos(3)/old_ax_width;
% end

figure(h_fig);%make the figure current
cb_axis = axes('Position',cb_pos);


for i = 1:3
    img(1:64,1,i) = linspace(c1(i),c2(i),64);
end
img(:,2,:) = img(:,1,:);
image(img,'Parent',cb_axis)

set(cb_axis,'XTick',[],'XTickLabel',[]);
set(cb_axis,'YAxisLocation','right','YDir','normal');
yl = get(cb_axis,'YLim');
ytick = linspace(yl(1),yl(2),11);
yticklabel = linspace(0,1,11);
set(cb_axis,'YTick',ytick, 'YTickLabel',yticklabel);

if ~isempty(ud)
    set(cb_axis,'YTickLabel',[]);
end

ud = get(ax,'UserData'); %stores colorbar handles;
ud = [ud, cb_axis];
set(ax,'UserData',ud);
if ~isempty(title)
    hl = ylabel(cb_axis,title);
    pos = get(hl,'Position');
    pos(1) = -1; %pos(1) - 1.1;
    set(hl,'Position',pos,'Interpreter','latex');
end


if ~isempty(varargin)
    str = [];
    for i = 1:length(varargin)
        if ischar(varargin{i})
            str = [str ',''' varargin{i} ''''];
        else
            str = [str ',' num2str(varargin{i}) ''];
        end
    end
    eval(['set(cb_axis' str ');']);
end


if strcmp(get(ax,'Type'),'axes')
    axes(ax);
end