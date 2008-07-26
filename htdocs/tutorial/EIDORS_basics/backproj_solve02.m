% $Id$

for idx=1:3
  if     idx==1; mdltype= 'b2c';
  elseif idx==2; mdltype= 'd2c2';
  elseif idx==3; mdltype= 'd2t2';
  end

   imdl= mk_common_model(mdltype,16);
   fmdl= imdl.fwd_model;
   params= aa_fwd_parameters( fmdl );

   homog= ones(size(fmdl.elems,1),1);
   img= eidors_obj('image','','elem_data',homog);
   img.fwd_model= fmdl;

   node_v= calc_all_node_voltages( img );
   elem_v= reshape(node_v( fmdl.elems',:),3,[],16);
   elem_v= squeeze(mean(elem_v,1));
   meas_v= params.N2E * node_v;

   sel= 4;
   ed= elem_v(:,sel);
   img.elem_data= ed;
   subplot(2,3,idx);
   show_fem(img);

   ej = zeros(size(ed));
   for i=1:16
     ej= ej + ( ed < meas_v(i,sel) );
   end
   img.elem_data= ej;
   subplot(2,3,idx+3);
   show_fem(img);

end

print -r125 -dpng backproj_solve02a.png;
