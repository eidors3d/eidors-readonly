load test_phess_solve_data.mat

%Chose a sep rad idx and plot
sepidx=1; radidx=1;

x0 = [-1/sqrt(2),-1/sqrt(2)];
xf = [1/sqrt(2),1/sqrt(2)];
img_pt_disc_line = ...
    line_plot_2d_img(img_pt_disc{sepidx,radidx},x0,xf,500);        
img_gn0_line = ...
    line_plot_2d_img(img_gn0{sepidx,radidx},x0,xf,500);       
img_eid_abs_line = ...
    line_plot_2d_img(img_eid_abs{sepidx,radidx},x0,xf,500);        


figure; hold on
plot(img_pt_disc_line,'r') %PT lbfp
plot(img_gn0_line,'b')
plot(img_eid_abs_line,'g')
