% Test cache deletes in the right order
% TODO: This needs to be implemented to test
%       that stuff gets deleted in priority order
function cache_limit_test

%   eidors_cache clear

    obj.type= 'data';
    for i=1:10
       obj.i=i;
       eidors_obj('set-cache',obj, 'cache_limit_test', i*ones(100));
       pause(.01);
    end

    disp('expect no output'); 
    eidors_cache show_objs

