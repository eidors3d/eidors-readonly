%$Id$
imdl= mk_common_model('c2c',16);

fmdl= imdl.fwd_model;
    ee= fmdl.elems';
    xx= reshape( fmdl.nodes(ee,1),3, []); 
    yy= reshape( fmdl.nodes(ee,2),3, []); 

homog= ones(size(fmdl.elems,1),1);
img= eidors_obj('image','','elem_data',homog);


for i=1:4
   if     i==1; stim= [0,1];
   elseif i==2; stim= [0,2];
   elseif i==3; stim= [0,4];
   elseif i==4; stim= [0,8];
   end

   fmdl.stimulation = mk_stim_patterns(16,1,stim,[0 1], {}, 1);

   img.fwd_model= fmdl;

   node_v= calc_all_node_voltages( img );

   % show voltages at measurement 2
   zz= reshape(     node_v(ee,4),3, []); 
   subplot(2,4,i); cla
   patch(xx,yy,zz,zz); view(0, 4); axis off
   subplot(2,4,i+4); cla
   patch(xx,yy,zz,zz); view(0,34); axis off
end

print_convert backproj_solve01a.png
