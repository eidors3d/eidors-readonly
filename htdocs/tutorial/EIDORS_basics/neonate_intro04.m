% Show plots
imgs = calc_slices(inv_solve(imdl, vh, vv));
axes('position',[0.32,0.6,0.63,0.25]);

taxis =  (0:size(imgs,3)-1)/13; % frame rate = 13
hold all
for i = 1:length(xposns);
    plot(taxis,ofs(i)+squeeze(imgs(yposns(i),xposns(i),:)),'LineWidth',2);
end
hold off
set(gca,'yticklabel',[]); xlim([0 16]);

print_convert neonate_intro04a.jpg
