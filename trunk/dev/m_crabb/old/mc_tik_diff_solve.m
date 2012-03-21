function [imgrec]=mc_tik_diff_solve(vh,vi,J,imghom,alphatik)
%Function that performs the most basic Tikhonov regularised solution for
%difference data in EIT. The output is an image strucutre with
%reconstructed conductivities

%Get the voltage difference vector, find total measurements
volthom=vh.meas; voltinhom=vi.meas; 
voltdiff=volthom-voltinhom;

%Now solver the linear system, using backslash (no regularisation)
condchange=tikhonov_reg(J,voltdiff,alphatik);

%Add the known homogeneous conductivity
condrecon=imghom.elem_data-condchange;

%Copy the image from homogeneous an store the elem data (condcutivities)
imgrec=imghom; imgrec.elem_data=condrecon;