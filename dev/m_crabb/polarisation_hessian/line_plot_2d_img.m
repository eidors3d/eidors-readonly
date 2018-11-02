function img_lp = line_plot_2d_img(img,x0,xf,n_points)
% LINE_PLOT_2D_IMG : Plot straight line through 2D image
% img_lp = line_plot_2d_img(img,x0,xf,n_points)
% img - 2D image
% x0/xf - start/end Cartesian coordinates (assumed lie on boundary): dim 2
% n_points is number of (equispaced) points for line plot through
% img_lp is image along 1D line: dim n_points

if(xf(1)==x0(1)) %vertical line plot
else
    m_ = (xf(2)-x0(2))/(xf(1)-x0(1)); c_ = xf(2)-m_*xf(1);
end

elem_list = img.fwd_model.elems;
node_list = img.fwd_model.nodes;

%Specific points along straight line
for i=1:n_points+1
    if(xf(1)==x0(1)) %vertical line plot just regular in y
        xi = x0(1);
        yi = x0(2) +(i-1)*(xf(2)-x0(2))/n_points;
        xycoords(i,:)=[xi,yi];    
    else
        xi = x0(1)+(i-1)*(xf(1)-x0(1))/n_points; %Regular in x
        yi = m_*xi + c_;
        xycoords(i,:)=[xi,yi];          
    end
end

%Find element for which each point belongs
t_lp = tsearchn(node_list,elem_list,xycoords);

%Plot the image
img_lp(1)=1; img_lp(n_points+1)=1; %assume background =1 
for i=2:n_points
    img_lp(i)=img.elem_data(t_lp(i));    
end

end
