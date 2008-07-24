%EIDORS 3D Matlab based package
%Version 2.0
%Matlab 6.1 is required for some of the plotting functions.
%
%-------------------------------------------------------------------------
%[1] function [E,D,Ela,pp] = fem_master_full(vtx,simp,mat,gnd_ind,elec,zc,perm_sym);
%
%Builds up the system matrix based on the complete electrode model. E is not 
%yet permuted. To permute E -> E(pp,pp) as in fwd_solver.
%
%--------------------------------------------------------------------------
%[2] function [Ef,D,Ela] = bld_master_full(vtx,simp,mat,elec,zc);
%
%System matrix assembling based on the complete electrode model.
%This function is called within fem_master_full.
%
%--------------------------------------------------------------------------
%[3] function [Ef,D,Ela] = bld_master(vtx,simp,mat_ref);
%
%Builds up the main compartment (GAP-SHUNT) of the system matrix for the complete
%electrode model. It is called within the function fem_master_full
%
%--------------------------------------------------------------------------
%[4] function [Er] = ref_master(E,vtx,gnd_ind,sch);
%
%Applys reference conditions to the system. Modifying the system matrix to preserve
%the uniqueness of the forward solution.
%
%--------------------------------------------------------------------------
%[5] function [I,Ib] = set_3d_currents(protocol,elec,vtx,gnd_ind,no_pl);
%
%This function sets current patterns in a system with (no_pl) planes of 
%equal number of electrodes according to "opposite" or "adjacent" protocols, 
%or their 3D similar.
%
%--------------------------------------------------------------------------
%[6] function [I,Ib] = set_multi_currents(protocol,elec,vtx,gnd_ind,no_pl);
%
%This functions applies opposite or adjacent current patterns to each of
%the planes of the system simultaneously. 
%
%--------------------------------------------------------------------------
%[7] function [V] = forward_solver(E,I,tol,pp,V);
%
%This function solves the forward problem using the Cholesky or LU method or 
%conjugate gradients. 
%
%--------------------------------------------------------------------------
%[8] function [voltageH,voltageV,indH,indV,df] = get_3d_meas(elec,vtx,V,Ib,no_pl);
%
%This function extracts multi-plane voltage measurements from a calculated
%3D nodal potential distribution V inside a tank with (no_pl) electrode
%planes. Each plane holds the same number of electrodes. Only the non-current
%caring electrodes at the time are involved in the measurements.
%
%--------------------------------------------------------------------------
%[9] function [voltage,ind,df] = get_multi_meas(protocol,elec,V,I,vtx,no_pl);
%
%The function can be used in the occasions where plane current patterns 
%(adjacent or polar) are adopted for systems with more planes, i.e. set by
%set_multi_currents function. Only non-current carrying electrodes are 
%involved in the measurements.
%
%--------------------------------------------------------------------------
%[10] [v_f] = m_3d_fields(vtx,el_no,m_ind,E,tol,gnd_ind,v_f);
%
%This function calculates the measurement fields using preconditioned
%conjugate gradients. These are used in the calculation of the Jacobian.
%
%--------------------------------------------------------------------------
%[11] function [fc] = figaro3d(srf,vtx,simp,fc,BB,h);
%
%This function plots the solution as a 3D object crossed at the plane z=h
%within the 3D boundaries of the volume.
%
%--------------------------------------------------------------------------
%[12] function [fc] = slicer_plot(h,BB,vtx,simp,fc);
%
%This function plots a 2D slice of the 3D solution vector BB at z=h.
%The function plots a smooth approximation to the calculated solution.
%
%--------------------------------------------------------------------------
%[13] function [fc] = slicer_plot_n(h,BB,vtx,simp,fc);
%
%This function plots a 2D slice of the 3D solution vector BB at z=h.
%
%--------------------------------------------------------------------------
%[14] function [CC] = solution_ext(BB,vtx,simp);
%
%Auxiliary function that extracts a secondary NODE-wise solution CC
%based on the calculated ELEMENT-wise solution BB. This function is
%called from slicer_plot function.
%
%--------------------------------------------------------------------------
%[15] function [J] = jacobian_3d(I,elec,vtx,simp,gnd_ind,mat_ref,zc,v_f,df,tol,perm_sym);
%
%This function calculates the Jacobian (sensitivity) matrix of the system wrt to conductivity.
%
%--------------------------------------------------------------------------
%[16] function [Jrr,Jri,Jir,Jii] = 
% jacobian_3d_comp(I,elec,vtx,simp,gnd_ind,mat_ref,no_pl,zc,v_f,df,tol,perm_sym);
%
%This function calculates the Jacobian matrices (wrt conductivity) for the complex EIT system.
%
%--------------------------------------------------------------------------
%[17] function [JTb] = adjoint_spin(vtx,simp,elec,x,gnd_ind,zc,I,no_pl,Vmes);
%
%The function calculates the product J'*b, i.e. Jacobian transpose times 
%a (measurements) vector b, using the adjoint sources formulation.
% 
%----------------------------------------------------------------------------
%[18] function [IntGrad] = integrofgrad(vtx,simp,mat_ref);
%
%function that calculates the integral of the gradients for first order
%tetrahedral elements. Required for the calculation of the Jacobian.
%
%--------------------------------------------------------------------------
%[19] function [elec_face,sels,cnts,VV] = set_electrodes(vtx,srf,elec_face,sels,cnts,VV);
%
%You need to call this function recursively to select boundary faces building up
%the ele_face. You will need to reshape this matrix appropriately to get the elec
%matrix, depending on how many faces there are in each electrode.
%
%--------------------------------------------------------------------------
%[20] function [mat,grp] = set_inho(srf,simp,vtx,mat_ref,val);
%
%Auxiliary functions used to set small local inhomogeneities
%in a volume (graphically).
%
%--------------------------------------------------------------------------
%[21] function [sel] = laserbeam(vtx,srf,cnts);
%
%Auxiliary plotting function
%Licenced
%--------------------------------------------------------------------------
%[22] function paint_electrodes(sel,srf,vtx);
%
%Auxiliary function which plots the electrodes red at the boundaries.
%
%--------------------------------------------------------------------------
%[23] function repaint_inho(mat,mat_ref,vtx,simp);
%
%Repaints the simulated inhomogeneity according to the reference
%distribution. (Increase -> Red, Decrease -> Blue) 
%
%-------------------------------------------------------------------------
%[24] function potplot(vtx,srf,V);
%
%Animates the forward solution in 3D.
%
%--------------------------------------------------------------------------
%[25] function [dd] = db23d(x1,y1,z1,x2,y2,z2);
%
%Auxiliary function that calculates the distance between two points in 3D.
%
%--------------------------------------------------------------------------
%[26] function [srf] = find_boundary(simp);
%
%Auxiliary function that calculates the boundary faces of a given 3D volume.
%Useful in electrode assignment.
%
%--------------------------------------------------------------------------
%[27] function [vtx_n,simp_n] = delfix(vtx,simp)
%
% Auxiliary function to remove the zero area faces
% produced by Matlab's delaunay triangulation
%
%---------------------------------------------------------------------------
%[28] function [vols] = check_vols(simp,vtx);
%
%Auxiliary function which calculates the volume of each tetrahedron in the mesh. 
%
%--------------------------------------------------------------------------
%[29] function [ta] = triarea3d(V);
%
%Function that calculates the area of a triangle in the 3D Cartesian space. 
%
%--------------------------------------------------------------------------
%[30] function [Reg] = iso_f_smooth(vtx,simp,deg,w);
%
%Calculates a first order discrete Gaussian smoothing operator
%
%--------------------------------------------------------------------------
%[31] function [Reg] = iso_s_smooth(simp,w);
%
%Calculates a second order discrete Gaussian smoothing operator
%
%--------------------------------------------------------------------------
%[32] function [solf,solp] = 
% inverse_solver(I,voltage,tol,mat_ref,vtx,simp,elec,no_pl,zc,perm_sym,gnd_ind,tfac,Reg,it);
%
%Calculates a Newton non-linear inverse solution.
%
%--------------------------------------------------------------------------
%[33] function demo_real( ... );
%
%Demonstrates an ERT (conductivity only) image reconstruction example.
%
%--------------------------------------------------------------------------
%[33] function demo_comp( ... );
%
%Demonstrates an EIT (conductivity and permittivity) image reconstruction example.
%
%--------------------------------------------------------------------------
%[34] function demo_complex(....);
%
%An alternative fully complex formulation for the complex EIT problem
%
%--------------------------------------------------------------------------
%-------------------------------------------------------------------------- 
