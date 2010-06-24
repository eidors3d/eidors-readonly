function mdl= create_gmsh_2d_circle(rad, n_elec)
% Create a 2D Circular FEM using GMSH
% mdl= CREATE_GMSH_2D_CIRCLE(rad, n_elec)
%
% mdl - EIDORS forward model
% rad - model radius

% License: GPL version 2 or version 3
% $Id$

p=[];

% target element size at each point
p=gmfl(p,'es = 0.009;');

theta= linspace(0,2*pi, n_elec+1); theta(end)=[];
for th=theta
  x=rad*sin(th);
  y=rad*cos(th);
  z=0;
  p=gmfn(p,'Point(%d) = {%f,%f,%f,es};', x,y,z);
end
  

fclose(p.fid);




% print to file - no number
function p=gmfl(p, str,  varargin)
   p= assign_p(p);
   fprintf(p.fid, [str,'\n'], varargin{:}); 

% print to file - add number
function p=gmfn(p, str,  varargin)
   p= assign_p(p);
   fprintf(p.fid, [str,'\n'], p.ptno, varargin{:}); 
   p.ptno = p.ptno + 1;

function p=assign_p(p);
   if isempty(p)
      % Unique name to about the second
      p.nam = sprintf('gmshmdl_%06d.geo', round(rem(now*1e5,1e6)));
      p.fid= fopen(p.nam,'wa');
      p.ptno= 0;
   end
