function ok= var_id_test
% Test that the mex file eidors_var_id works 
% $Id: var_id_test.m,v 1.1 2005-10-25 13:40:07 aadler Exp $

ok=1;

%
% Test 2:
%   test for random sorting
%
n_var= 100;
for i=1:n_var
   str{i} = char( 'a' + floor(26*rand(1,20)) );
end

% Randomly order them, and assign to a variable
var_id= '';
for iter=1:100
   [jnk, idx] = sort( rand(1,n_var) );
   
   vv= struct([]);
   for i= 1:n_var; %idx
      eval(['vv(1).', str{i}, '=i;']);
   end

   evi= eidors_var_id( vv );
   if isempty(var_id)
      var_id= evi;
   elseif ~strcmp( var_id, evi );
         ok=0;
keyboard
   end
end



