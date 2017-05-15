ne=length(img(1).fwd_model.electrode);
clf; h=show_fem(img(1).fwd_model,[0 1]); h.EdgeColor=[1 1 1];
axis off; set(gcf,'Color',[1 1 1]);
t=get(gca,'Children');
for i = [1 3 5 7 9 11];
   t(i).FontSize = 40; t(i).Color = [0 0 0];
end
for i=[1:ne];
   n = img(1).fwd_model.electrode(i).nodes;
   xy(i,:) = mean(img(1).fwd_model.nodes(n,:));
end
for i=[1:ne]; for j=[1:ne];
   line(xy([i j],1),xy([i j],2));
end; end
print_convert('model_reduction02.png');
