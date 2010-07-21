function img= bs_nonlinearGN( inv_model, data )
% NONLINEARGN Non-Linear EIT Inverse
% img        => output image (or vector of images)
% inv_model  => inverse model struct
% data       => measurement data
% $Id: nonlinearGN.m 1535 2008-07-26 15:36:27Z aadler $
% Modified by Bartosz Sawicki

fwd_model= inv_model.fwd_model;

img= eidors_obj('image','Solved by bs_nonlinearGN');
img.fwd_model= fwd_model;
sol= inv_model.nonlinearGN.init_backgnd * ...
     ones(size(fwd_model.coarse2fine,2),1);

RtR = calc_RtR_prior( inv_model );
hp2 = calc_hyperparameter( inv_model )^2;

factor= 0; norm_d_data= inf;

for iter= 1:inv_model.parameters.max_iterations
   img.elem_data= sol;
   simdata= fwd_solve( img );

   d_data= data - simdata.meas;
   prev_norm_d_data= norm_d_data; norm_d_data= norm(d_data); 
   eidors_msg('bs_nonlinearGN: iter=%d diff=%f factor=%f', ...
           iter, norm_d_data, factor, 2);

   if prev_norm_d_data - norm_d_data < inv_model.parameters.term_tolerance
      eidors_msg('   BREAK (Requested tolerance achieved)', 2);
      break;
   end
   
   J = calc_jacobian( fwd_model, img);
   delta_sol = (J.'*J + hp2*RtR)\ (J.' * d_data);
     
   factor = linesearch(img, data, sol, delta_sol, norm_d_data);

   sol = sol + factor*delta_sol;
end

   img.elem_data= sol;

% Test several different factors to see which minimizes best
function factor= linesearch(img, data, sol, delta_sol, norm_d_data)
   facts= linspace(0,2,20);  
   norms= norm_d_data;
   for f= 1:length(facts);
      img.elem_data= sol + facts(f)*delta_sol;
      simdata= fwd_solve( img );
      norms(f)= norm(data - simdata.meas);
   end
   ff= find(norms==min(norms));
   factor= facts(ff(end));
