imdl = mk_common_model('c2C',16); img = mk_image(imdl.fwd_model,1);
vv=fwd_solve(img);      v0e=vv.meas;
JJ=calc_jacobian(img);  J04=JJ(4,:)';

%High-order EIDORS solver %Change default eidors solvers
img.fwd_model.solve = @fwd_solve_higher_order;
img.fwd_model.system_mat = @system_mat_higher_order;
img.fwd_model.jacobian = @jacobian_adjoint_higher_order;

vve=[]; JJ4=[];
for i= 1:3; switch i;
   case 1; img.fwd_model.approx_type = 'tri3'; % linear
   case 2; img.fwd_model.approx_type = 'tri6'; % quadratic
   case 3; img.fwd_model.approx_type = 'tri10'; % cubic;
   end %switch
   vv=fwd_solve(img);      vve(:,i)=vv.meas;
   JJ=calc_jacobian(img);  JJ4(:,i)=JJ(4,:)';
end

subplot(311);
plot([v0e,vve,(v0e*[1,1,1]-vve)*10]);
legend('Default','linear','quadratic','cubic','(1-0)x10','(2-0)x10','(3-0)x10');
xlim([1,100]);

imgJJ=img; imgJJ.elem_data = JJ4;
imgJJ.show_slices.img_cols = 3;

subplot(312); show_slices(imgJJ); eidors_colourbar(imgJJ);

imgJJ.elem_data = JJ4 - J04*[1,1,1];
subplot(313); show_slices(imgJJ); eidors_colourbar(imgJJ);

print_convert forward_solvers_2d_high_order05a.png
