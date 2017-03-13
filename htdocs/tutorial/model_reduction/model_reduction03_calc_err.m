function err = model_reduction03_calc_err(RR, mm, tol)
   err = zeros(size(mm));
   for ii = 1:length(mm);
      mi = mm(ii);
      err_mm = ((model_reduction03_stdres(RR*mi,tol)-RR*mi)-1)/mi;
      err(ii) = norm(err_mm(:));
   end
