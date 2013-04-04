% $Id$
xl=-3; xr= 3; yb=-15; yt= 15;
np= 35;
[x,y] = meshgrid( linspace(xr,xl,np), linspace(yb,yt,61) );
vtx= [y(:),x(:)];
for i=1:np
   elec_nodes{i   }= [y(1,i);x(1,i)]';
   elec_nodes{i+np}= [y(end,i);x(end,i)]';
end
z_contact= 1e5;
fmdl= mk_fmdl_from_nodes( vtx, elec_nodes, z_contact, 'sq_m1');
fmdl.stimulation(1).stimulation='Amp';
fmdl.stimulation(1).stim_pattern=[ones(np,1);-ones(np,1)];
fmdl.stimulation(1).meas_pattern=zeros(1,2*np); % don't care
fmdl2=fmdl;


%Pick some elements in a rectangular box
idx=[]; a=0; n_elems=size(fmdl2.elems,1);

%Define a rectangular region
for i=1:n_elems
    enodes = fmdl2.elems(i,:); 
    nd = mean(fmdl2.nodes(enodes,:)); %centroid
    %Select upper rectanular box in x and y
    if(nd(1) < 14 && nd(1) >10 && nd(2) <1.55 && nd(2) >-1.55 )
        a=a+1;
        idx=[idx,i];
    end                        
end

%Isotropic material in box
img = mk_image( fmdl, 100);
img.elem_data(idx)=1;

% Solve forward problem for isotropic
img.fwd_solve.get_all_meas = 1;
vh = fwd_solve(img);
img.fwd_model.mdl_slice_mapper.npx = 200;
img.fwd_model.mdl_slice_mapper.npy = 100;
subplot(221);
q = show_current(img,vh.volt);
hh=show_fem(img);
set(hh,'EdgeColor',[1,1,1]*.75);
hold on;
sy = linspace(-3,3,25); sx =  15+0*sy;
hh=streamline(q.xp,q.yp, q.xc, q.yc, sx,sy); set(hh,'Linewidth',2);
hold off;

%Equipotentials on isotropic
imge = rmfield(img, 'elem_data');
imge.node_data = vh.volt(:,1);
subplot(222); show_fem(imge);

%Anisotropic material in box
subplot(223); show_fem(img)
%Assign the anisotropic tensor
for i=1:n_elems
    elem_data(i,1,1,1) =100;
    elem_data(i,1,1,2) =0; elem_data(i,1,2,1) = elem_data(i,1,1,2);
    elem_data(i,1,2,2) =100;       
end
%Anisotropic part for ractangle
for i=idx
    elem_data(i,1,1,1) =1;
    elem_data(i,1,1,2) =0.99; elem_data(i,1,2,1) = elem_data(i,1,1,2);
    elem_data(i,1,2,2) =1;       
end

%Solve forward problem for isotropic
fmdl.system_mat = @system_mat_higher_order_anisotropy;
fmdl.approx_type    = 'tri3'; % linear
img2 = mk_image_anisotropy(fmdl,elem_data);
img2.fwd_solve.get_all_meas = 1; %Internal voltage
v2 = fwd_solve(img2); 
img2.fwd_model.mdl_slice_mapper.npx = 200;
img2.fwd_model.mdl_slice_mapper.npy = 100;
subplot(223);
q = show_current(img2,v2.volt);
hh=show_fem(img);
set(hh,'EdgeColor',[1,1,1]*.75);
hold on;
sy = linspace(-3,3,25); sx =  15+0*sy;
hh=streamline(q.xp,q.yp, q.xc, q.yc, sx,sy); set(hh,'Linewidth',2);
hold off;

%Equipotentials on anisotropic
imgp = rmfield(img2, 'elem_data');
imgp.node_data = v2.volt(:,1);
subplot(224); show_fem(imgp);