function img= tutorial151_nonlinearGN( inv_model, data )
% TUTORIAL151_NONLINEARGN Non-Linear EIT Inverse
% img        => output image (or vector of images)
% inv_model  => inverse model struct
% data       => measurement data
% $Id: tutorial151_nonlinearGN.m,v 1.1 2006-11-17 03:23:29 aadler Exp $

fwd_model= inv_model.fwd_model;

img= eidors_obj('image','Solved by tutorial151_nonlinearGN');
img.fwd_model= fwd_model;
sol= inv_model.tutorial151_nonlinearGN.init_backgnd * ...
     ones(size(fwd_model.elems,1),1);

    RtR = calc_RtR_prior( inv_model );
    hp2  = calc_hyperparameter( inv_model )^2;

for iter= 1:inv_model.parameters.max_iterations
   img.elem_data= sol;
keyboard
   simdata= fwd_solve( img );

   data_diff= data - simdata.meas;
   eidors_msg('tutorial151_nonlinearGN: iter=%d err=%f', ...
           iter,norm(data_diff), 3);
   if norm(data_diff) < inv_model.parameters.term_tolerance
      break;
   end
   

   J = calc_jacobian( fwd_model, img);
   delta_sol = (J.'*J + hp2*RtR)\ (J.' * data_diff);
   sol = sol + delta_sol;
end

   img.elem_data= sol;



