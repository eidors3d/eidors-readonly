function [index_incyl,mat]=elems_in_cylinder(center_simp,origin,radius,mat,value);
% INCYL: calculate elements in infinite cylinder
% [index_incyl,mat]=incyl(center_simp,origin,radius,mat,value);
%
% See also: mk_c2f_circ_mapping (this is the recommended function)
 
% (C) 2005 David Stephenson. Licenced under GPL
% $Id$

index_incyl=[];
a=origin(1);
b=origin(2);
r=radius;

for i=1:size(center_simp,1);
    
	x=center_simp(i,1);
	y=center_simp(i,2);	 
	r_squared_calc=((x-a)^2)+((y-b)^2);
	
    if r_squared_calc <= r^2
        index_incyl(i,1)=1;
	else
        index_incyl(i,1)=0;
	end    
       
    i=i+1;
    
end

for i=1:size(center_simp);
        
	if index_incyl(i,1)==1;
        
        mat(i,1)=value;
        
        i=i+1;
        
	end
    
end
