function ok= var_id_test
% Test that the mex file eidors_var_id works 

% (C) 2005 Andy Adler. License: GPL version 2 or version 3
% $Id: var_id_test.m,v 1.8 2007-08-29 09:07:16 aadler Exp $

ok=1;

%
% Test 1:
%  Test variable types
%
vv1.a=1; vv1.b='asfd'; vv1.c(1)=1; vv1.c(2)=2; vv1.s= @sin;
vv2.a=1; vv2.b='asfd'; vv2.c(1)=1; vv2.c(2)=2; vv2.s= @sin;
if ~strcmp( eidors_var_id(vv1), eidors_var_id(vv2) )
   warning('var_id_test: 1');
   ok=0;
end

vv2.a=1; vv2.b='asfd'; vv2.c(1)=1; vv2.c(2)=2; vv2.s= @sin;
vv2.c(2)=3;
if strcmp( eidors_var_id(vv1), eidors_var_id(vv2) )
   warning('var_id_test: 2'); ok=0;
end

vv2.a=1; vv2.b='asfd'; vv2.c(1)=1; vv2.c(2)=2; vv2.s= @sin;
vv2.b='asdf';
if strcmp( eidors_var_id(vv1), eidors_var_id(vv2) )
   warning('var_id_test: 3'); ok=0;
end

vv2.a=1; vv2.b='asfd'; vv2.c(1)=1; vv2.c(2)=2; vv2.s= @sin;
vv2.s= @cos;
if strcmp( eidors_var_id(vv1), eidors_var_id(vv2) )
   warning('var_id_test: 4'); ok=0;
end

vv2.a=1; vv2.b='asfd'; vv2.c(1)=1; vv2.c(2)=2; vv2.s= @sin;
vv2.a=1 + eps;
if strcmp( eidors_var_id(vv1), eidors_var_id(vv2) )
   warning('var_id_test: 5'); ok=0;
end


%
% Test 2:
%   test for random sorting
%
n_var= 10;
for i=1:n_var
   str{i} = char( 'a' + floor(26*rand(1,20)) );
end

% Randomly order them, and assign to a variable
var_id= '';
for iter=1:2
   [jnk, idx] = sort( rand(1,n_var) );
   
   vv= struct([]);
   for i= 1:n_var; %idx
      eval(['vv(1).', str{i}, '=i;']);
   end

   evi= eidors_var_id( vv );
   if isempty(var_id)
      var_id= evi;
   elseif ~strcmp( var_id, evi );
      warning('var_id_test: 6'); ok=0;
   end
end



