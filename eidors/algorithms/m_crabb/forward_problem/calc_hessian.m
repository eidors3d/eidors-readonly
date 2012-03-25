function H = calc_hessian(fwd_model,img)
%Find the Hessian associated with an image (and forward model)
%Second derivative of discretization method
%This is just a cached version of code - mc_calc_hessian

cache_obj= {fwd_model,img};
H = eidors_obj('get-cache', cache_obj, 'hessian');
if ~isempty(H)
            eidors_msg('hessian: using cached value', 3);        
        return
end

[H]= feval(fwd_model.hessian, fwd_model, img);

eidors_obj('set-cache', cache_obj, 'hessian',H);
eidors_msg('hessian: setting cached value', 3);

end