function [bound,elem,nodes]=fem_modify(fwd_model)
% FEM_MODIFY:  Modify the FEM for high order FEM called as
%    [bound,elem] = mc_fem_modify( fwd_model )
% It will call :fwd_model.modify=@mc_fem_modify
% where we need to make sure fwd_model.mc_type='tri3', 'tet4' etc. 
%
% fwd_model : is a fwd_model structure
% bound : is the boundary nodes numbers (from boundary (bound(ii).nodes))
% elem : is the element nodes numbers (from elems elem(ii).nodes))

cache_obj= {fwd_model};
bound = eidors_obj('get-cache', cache_obj, 'bound');
elem = eidors_obj('get-cache', cache_obj, 'elem');
nodes=eidors_obj('get-cache',cache_obj,'nodes');
if ~isempty(bound)
    if ~isempty(elem)
        if ~isempty(nodes)
            eidors_msg('bound: using cached value', 3);        
            eidors_msg('elem: using cached value',3);
            eidors_msg('nodes: using cached value',3);
        return
        end
    end
end

[bound,elem,nodes]= feval(fwd_model.fem_modify, fwd_model);

eidors_obj('set-cache', cache_obj, 'bound', bound);
eidors_obj('set-cache', cache_obj, 'elem',elem);
eidors_obj('set-cache', cache_obj, 'nodes',nodes);
eidors_msg('fem_modify: setting cached value', 3);