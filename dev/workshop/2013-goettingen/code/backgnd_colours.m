bkgnds = [ ...
    50,50,15;
    40,40,20;
    30,30,20;
    40,40,30;
    40,40,35;
    30,30,30;
    20,20,20;
    15,50,50;
    20,40,40;
    20,30,30;
    30,40,50;
    
    ]/100; 

imgs = eidors_obj('image','','elem_data',-10:1:10);
imgs.calc_colours.cb_shrink_move = [2,.8,.08];
for i=1:size(bkgnds,1);
    for j=1:3
    img(1,1,j) = bkgnds(i,j);
    end


    subplot(4,3,i);
    image(img)
    eidors_colourbar(imgs);
    
end
