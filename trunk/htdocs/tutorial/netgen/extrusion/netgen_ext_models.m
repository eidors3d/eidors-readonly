function netgen_ext_models( numbers )

if nargin==0; numbers = 1:8; end
for i=numbers(:)';
   do_sim(i);
end

function do_sim( number )
switch number

   case 1;
xy= [ -0.89 -0.74 -0.21  0.31  0.79  0.96  0.67  0.05 -0.36 -0.97;
       0.14  0.51  0.35  0.50  0.27 -0.23 -0.86 -0.69 -0.85 -0.46]';
[fmdl, mat_idx] = ng_mk_extruded_model({2,xy,1},[9,0,1],[0.02]);

   case 2;
xy= [ -0.89 -0.74 -0.21  0.31  0.79  0.96  0.67  0.05 -0.36 -0.97;
       0.14  0.51  0.35  0.50  0.27 -0.23 -0.86 -0.69 -0.85 -0.46]';
[fmdl, mat_idx] = ng_mk_extruded_model({2,xy,1,0.1},[5,0,1],[0.10]);

   case 3;
xy= [ -0.89 -0.74 -0.21  0.31  0.79  0.96  0.67  0.05 -0.36 -0.97;
       0.14  0.51  0.35  0.50  0.27 -0.23 -0.86 -0.69 -0.85 -0.46]';
[fmdl, mat_idx] = ng_mk_extruded_model({2,xy,[3,50]},[9,0,1],[0.05]);

   case 4;
xy= [ -0.89 -0.74 -0.21  0.31  0.79  0.96  0.67  0.05 -0.36 -0.97;
       0.14  0.51  0.35  0.50  0.27 -0.23 -0.86 -0.69 -0.85 -0.46]';
[fmdl, mat_idx] = ng_mk_extruded_model({1,xy,[3,50]},[9,0,1],[0.05]);


   case 9;
xy= [ -0.89 -0.74 -0.21  0.31  0.79  0.96  0.67  0.05 -0.36 -0.97;
       0.14  0.51  0.35  0.50  0.27 -0.23 -0.86 -0.69 -0.85 -0.46]';
[fmdl, mat_idx] = ng_mk_extruded_model({2,{a,0.3*a},1},[16,0,1],[0.01]);

end

show_fem(fmdl); view(170,20);
print_convert( ...
   sprintf('netgen_ext_models%02d_1.png',number), '-density 75');
view(0,90);
print_convert( ...
   sprintf('netgen_ext_models%02d_2.png',number), '-density 75');
