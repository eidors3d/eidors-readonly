function [H_obj,GN_only_obj,H_only_obj] = calc_hessian_diag_obj(fwd_model,img,elem_list,delta_d)
%Find the Hessian of objective function associated
%Second derivative of discretization method
%3rd input elem_list - indices of elementn we want Hessian for
%M Crabb - 29.06.2012
%
%H_obj = J'*J(:,:) + H(k,:,:)*dv(k)
%      = GN_only_obj + H_only_obj

%{
cache_obj= {fwd_model,img};
H_obj = eidors_obj('get-cache', cache_obj, 'hessian');
if ~isempty(H_obj)
            eidors_msg('hessian: using cached value', 3);        
        return
end
%}
J = calc_jacobian(fwd_model,img);
Je = J(:,elem_list);
[H]= calc_hessian_diag(fwd_model,img,elem_list);

GN_only_obj = Je'*Je;

H_only_obj = zeros(size(GN_only_obj));
for k=1:length(delta_d)
H_only_obj = H_only_obj + reshape(H(k,:,:)*delta_d(k),size(H,2),size(H,3));
end

H_obj = GN_only_obj + H_only_obj;

%{
eidors_obj('set-cache', cache_obj, 'hessian',H);
eidors_msg('hessian: setting cached value', 3);
%}

end