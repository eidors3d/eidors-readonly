function dg_show_gallery_models(image1,image2,title1,title2,title3)

image= image1;
% Convert to resistivity
image1.elem_data= 1./image1.elem_data;
image2.elem_data= 1./image2.elem_data;
image.elem_data= image1.elem_data - image2.elem_data;
m1= mean(image1.elem_data);
m2= mean(image2.elem_data);
m= round(mean([m1 m2])/10)*10;
% Draw both models
limit= 20;
calc_colours('npoints',128);
calc_colours('ref_level',m);
figure;
    show_slices(image1,[inf,inf,0],limit,m)
    calc_colours(image1,limit,1);
    title({title1},'FontSize',24);

figure;
    nodes= image1.fwd_model.nodes;
    x=nodes(:,1); y=nodes(:,2); z=nodes(:,3); nodes= [x z y];    
    image1.fwd_model.nodes= nodes;
    show_fem(image1,[1,0,limit]);
    calc_colours(image1,limit,1);
    title({title1},'FontSize',24);
    camproj('perspective');
    axis square;

figure;
    show_slices(image2,[inf,inf,0],limit,m);
    calc_colours(image1,limit,1);
    title({title2},'FontSize',24);
    
figure;
    nodes= image2.fwd_model.nodes;
    x=nodes(:,1); y=nodes(:,2); z=nodes(:,3); nodes= [x z y];    
    image2.fwd_model.nodes= nodes;
    show_fem(image2,[1,0,limit]);
    calc_colours(image1,limit,1);
    title({title2},'FontSize',24);
    camproj('perspective');
	axis square;
    
% Draw model difference
calc_colours('ref_level',0);
limit= 10;
figure;
    show_slices(image,[inf,inf,0],limit,0);
    calc_colours(image,limit,1);
    title({title3},'FontSize',24);

figure;
    nodes= image.fwd_model.nodes;
    x=nodes(:,1); y=nodes(:,2); z=nodes(:,3); nodes= [x z y];    
    image.fwd_model.nodes= nodes;
    show_fem(image,[1,0,limit]);
    calc_colours(image,limit,1);
    title({title3},'FontSize',24);
    camproj('perspective');
    axis square;
    
end