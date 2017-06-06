nmesh = 3; %greater or equal 4 for sensible rate estimates (depends on computer architec)

%Build a square CEM model
%Nodal coordinates and electrode nodes 
x=linspace(0,pi,6); 
y=linspace(0,pi,6);
[X,Y]=meshgrid(x,y); vtx=[X(:) Y(:)];
elec_nodes={[0 0],[0,pi]};
mdl=mk_fmdl_from_nodes(vtx,elec_nodes,0.02,'name');

%Electrode node indices and nodes on boundary not on elecs
mdl.electrode(1).nodes = [7,13];
mdl.electrode(2).nodes = [19,25];
bound_nodes_not_elecs = [1,2,3,4,5,6,12,18,24,30,31,32,33,34,35,36];
    
%Make image and add some fields
img=mk_image(mdl,1); %figure; show_fem(img);
stim = mk_stim_patterns(2,1,'{ad}','{ad}',{'meas_current'}, 1);
mdl.stimulation=stim; 
img.fwd_model.bound_nodes_not_elecs=bound_nodes_not_elecs;   
mdl.bound_nodes_not_elecs=bound_nodes_not_elecs;
img.fwd_model.stimulation=stim; img.fwd_solve.get_all_meas=1;  
mdl_init=mdl;


%Copy model, assign unit conductivity and plot
mdl_h1=mdl; img_h1=mk_image(mdl_h1,1);
figure; 
subplot(121); show_fem(img_h1); hold on;
plot(img_h1.fwd_model.nodes(bound_nodes_not_elecs,1),...
img_h1.fwd_model.nodes(bound_nodes_not_elecs,2),'r*');
xlab = xlabel('x'); set(xlab,'FontSize',14);
ylab = ylabel('y'); set(ylab,'FontSize',14);      

%Refine model
mdl_h2=h_refine_square_domain(mdl); img_h2=mk_image(mdl_h2,1);
subplot(122); show_fem(img_h2); hold on;
plot(img_h2.fwd_model.nodes(bound_nodes_not_elecs,1),...
img_h2.fwd_model.nodes(bound_nodes_not_elecs,2),'r*');
xlab = xlabel('x'); set(xlab,'FontSize',14);
ylab = ylabel('y'); set(ylab,'FontSize',14);      

print_convert error_rates_contact_impedance01a.png
