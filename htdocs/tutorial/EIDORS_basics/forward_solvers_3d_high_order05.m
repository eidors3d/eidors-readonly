imdl = mk_common_model('b3cr',16); img = mk_image(imdl.fwd_model,1);
vv=fwd_solve(img);      v0e=vv.meas;
JJ=calc_jacobian(img);  J04=JJ(4,:)';

%High-order EIDORS solver %Change default eidors solvers
img.fwd_model.solve = @fwd_solve_higher_order;
img.fwd_model.system_mat = @system_mat_higher_order;
img.fwd_model.jacobian = @jacobian_adjoint_higher_order;

vve=[]; JJ4=[];
for i= 1:2; switch i;
   case 1; img.fwd_model.approx_type = 'tet4'; % linear
   case 2; img.fwd_model.approx_type = 'tet10'; % quadratic
   end %switch
   vv=fwd_solve(img);      vve(:,i)=vv.meas;
   JJ=calc_jacobian(img);  JJ4(:,i)=JJ(4,:)';
end

subplot(311);
plot([v0e,vve,(v0e*[1,1]-vve)*10]);
legend('Default','linear','quadratic','(1-0)x10','(2-0)x10');
xlim([1,100]);

imgJJ=img; imgJJ.elem_data = JJ4;
imgJJ.show_slices.img_cols = 2;

level = [inf,inf,0.3];
subplot(312); show_slices(imgJJ,level); eidors_colourbar(imgJJ);

imgJJ.elem_data = JJ4 - J04*[1,1];
subplot(313); show_slices(imgJJ,level); eidors_colourbar(imgJJ);

print_convert forward_solvers_3d_high_order05a.png
