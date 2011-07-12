%% calculate the GREIT parameters
levels =[inf,inf,Zpos];
show_slices(imgr, levels);
imgr.calc_colours.npoints = 128;
imgr.calc_slices.levels=levels;
params = eval_GREIT_fig_merit(imgr, xyzr);

p_names = {'AR','PE','RES','SD','RNG'};
for i=1:5; subplot(5,1,i);
    plot(params(i,:)); ylabel(p_names{i});
end
print_convert test_params04a.png '-density 100'
