function [center_simp]=center_of_simps(fwd_model, extraparam)
% CENTER_OF_SIMPS: Calculates the Center of Mass of the Simplicies.
%
% Usage type #1
% [center_simp]=center_of_simps(fwd_model);
% Usage type #2
% [center_simp]=center_of_simps(simp,vtx);
%
% center_simp = The center of mass of simps 
%               [Nelems x 3] (x y z co-ordinates)
%
% (C) 2005 David Stephenson. Licensed under GPL Version 2
% $Id: center_of_simps.m,v 1.4 2007-08-29 09:13:52 aadler Exp $

if nargin==1
    simp= fwd_model.elems;
    vtx = fwd_model.nodes;
else
    simp= fwd_model;
    vtx= extraparam;
end

    x_sum=vtx(simp(:,1),1)+vtx(simp(:,2),1)+vtx(simp(:,3),1)+vtx(simp(:,4),1);
    y_sum=vtx(simp(:,1),2)+vtx(simp(:,2),2)+vtx(simp(:,3),2)+vtx(simp(:,4),2);
    z_sum=vtx(simp(:,1),3)+vtx(simp(:,2),3)+vtx(simp(:,3),3)+vtx(simp(:,4),3);
    
    x_center=x_sum/4;
    y_center=y_sum/4;
    z_center=z_sum/4;
    
    center_simp(:,1)=x_center;
    center_simp(:,2)=y_center;
    center_simp(:,3)=z_center;
    


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is part of the EIDORS suite.
% Copyright (c) D.R Stephenson 2004
% Copying permitted under terms of GNU GPL
% See enclosed file gpl.html for details.
% EIDORS 3D version XXX
% MATLAB Version 6.5.0.180913a (R13)
% MATLAB License Number: 1560
% Operating System: Microsoft Windows XP Version 5.1 (Build 2600: Service Pack 1)
% Java VM Version: Java 1.3.1_01 with Sun Microsystems Inc. Java HotSpot(TM) Client VM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
