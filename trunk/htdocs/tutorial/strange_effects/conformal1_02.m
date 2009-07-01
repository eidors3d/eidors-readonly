vv = fwd_solve(img);
node_v= vv.volt;

ee= img.fwd_model.elems';
xx= reshape( img.fwd_model.nodes(ee,1),3, []); 
yy= reshape( img.fwd_model.nodes(ee,2),3, []); 

subplot(222); cla;
   % show voltages at measurement 2
   zz= reshape(     node_v(ee,1),3, []); 
   patch(xx,yy,zz,zz); view(0, 90); axis image
%  subplot(2,4,i+4); cla
%  patch(xx,yy,zz,zz); view(0,34); axis off

vv = fwd_solve(img2);
node_v= vv.volt;

ee= img2.fwd_model.elems';
xx= reshape( img2.fwd_model.nodes(ee,1),3, []); 
yy= reshape( img2.fwd_model.nodes(ee,2),3, []); 

subplot(224); cla;
   % show voltages at measurement 2
   zz= reshape(     node_v(ee,1),3, []); 
   patch(xx,yy,zz,zz); view(0, 90); axis image
