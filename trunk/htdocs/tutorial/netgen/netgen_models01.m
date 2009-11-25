% Simple 3D cylinder. Radius = 1. No electrodes
    fmdl= ng_mk_cyl_models(3,[0],[]); 
% Simple 2D cylinder. Radius = 2. Set minsize to refine
    fmdl= ng_mk_cyl_models([0,2,.2],[0],[]); 
% 3D cylinder. Radius = 1. 2 planes of 8 elecs with radius 0.1
    fmdl= ng_mk_cyl_models(3,[8,1,2],[0.1]); 
% 3D cylinder. Radius = 1. 6 circ elecs with elec refinement
    fmdl= ng_mk_cyl_models(3,[7,1],[0.2,0,0.05]); 
% 3D cylinder. Radius = 1. 5 rect elecs with no refinement
    fmdl= ng_mk_cyl_models(3,[7,1],[0.2,0.3]); 
% 2D cylinder. Radius = 1. 11 rect elecs with refinement
    fmdl= ng_mk_cyl_models(0,[11],[0.2,0,0.05]); 
% 2D cylinder. Radius = 1.5. Refined(0.1). 11 elecs with refinement
    fmdl= ng_mk_cyl_models([0,1,0.1],[11],[0.2,0,0.02]); 
% 2D cylinder. elecs at 0, 90 and 120 degrees
    fmdl= ng_mk_cyl_models(0,[0;90;120],[0.2,0,0.03]); 
% 2D cylinder. elecs at 0 (large,refined) and 120 (small) degrees
    fmdl= ng_mk_cyl_models(0,[0;90],[0.4,0,0.01;0.1,0,0.1]); 
% 3D cylinder. elecs at 0, 30, 60, 90 in planes
    fmdl= ng_mk_cyl_models(3,[0,0.5;30,1;60,1.5;90,2.0],[0.2,0,0.1]); 
% 3D cylinder. Various elecs at 0, 30, 60, 90 in planes
    el_pos = [0,0.5;30,1;60,1.5;90,2.0];
    el_sz  = [0.2,0,0.1;0.1,0,0.05;0.2,0.2,0.02;0.2,0.4,0.5];
    fmdl= ng_mk_cyl_models(3,el_pos,el_sz); 
% Simple 3D cylinder with a ball
    extra={'ball','solid ball = sphere(0.5,0.5,2;0.4);'}
    [fmdl,mat_idx]= ng_mk_cyl_models(3,[0],[],extra); 
    img= eidors_obj('image','ball'); img.fwd_model= fmdl;
    img.elem_data(mat_idx{1}) = 1; img.elem_data(mat_idx{2}) = 2;
% 3D cylinder with 8 electrodes and cube
    extra={'cube','solid cube = orthobrick(0.5,0.5,0.5;0,0,1.5);'}
    [fmdl,mat_idx]= ng_mk_cyl_models(2,[8,0.5,1.5],[0.1],extra); 
% 3D cylinder with inner cylinder
    extra={'ball','solid ball = cylinder(0.2,0.2,0;0.2,0.2,1;0.2) and orthobrick(-1,-1,1;1,1,2) -maxh=0.05;'}
    [fmdl,mat_idx]= ng_mk_cyl_models(3,[0],[],extra); 
% 2D cylinder with 8 electrodes and hole
    extra={'ball','solid ball = sphere(0.2,0.2,0;0.2) -maxh=0.05;'}
    fmdl= ng_mk_cyl_models(0,[8],[0.1,0,0.05],extra); 
% 2D cylinder with 9 electrodes and inner cylinder
    extra={'ball','solid ball = cylinder(0.2,0.2,0;0.2,0.2,1;0.2) and orthobrick(-1,-1,0;1,1,0.05) -maxh=0.03;'}
    fmdl= ng_mk_cyl_models(3,[9],[0.2,0,0.05],extra); 
    img= eidors_obj('image','ball'); img.fwd_model= fmdl;
    ctr = interp_mesh(fmdl); ctr=(ctr(:,1)-0.2).^2 + (ctr(:,2)-0.2).^2;
    img.elem_data = 1 + 0.1*(ctr<0.2^2);

