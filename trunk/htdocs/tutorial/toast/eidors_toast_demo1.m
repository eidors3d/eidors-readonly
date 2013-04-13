% This functions needs the functions from the matlab 
%  version of toast to be installed. You can get them from 
%    http://web4.cs.ucl.ac.uk/research/vis/toast/
%
imdl=mk_common_model('a2c2',8);
imdl.fwd_model.stimulation = mk_stim_patterns(8,1,'{mono}','{ad}',{}',0.1);
fmdl=imdl.fwd_model;

nel = size(fmdl.elems,1);
vtx= fmdl.nodes; tri = fmdl.elems;
for i= 1:nel
    n1= vtx(tri(i,1),:);
    n2= vtx(tri(i,2),:);
    n3= vtx(tri(i,3),:);
    if det([n2-n1;n3-n1]) < 0;
       tri(i,:) = tri(i,[2 1 3]);
    end
end
subplot(221);show_fem(fmdl);
hMesh = toastMakeMesh (vtx,tri,15*ones(size(fmdl.elems,1)));

subplot(222); toastShowMesh(hMesh)
write_toast_qm('a2c2_8.qm',fmdl);
toastReadQM(hMesh,'a2c2_8.qm');


qvec = toastQvec(hMesh,'Neumann','Point');

% DOT Problem is -Grad k dot Grad phi + mua phi = 0
% EIT Problem is -Grad s dot Grad phi = 0
% k = 1/3/(mua + mus).
% To match EIT, set mua = 0
% s = k = 1/3/mus. mus = s/3
% ref is refractive index -> set to 1 for EIT
% BUT, we need to multiply k by 0.3 (mm/ps) (speed of light)
mua=zeros(nel,1);
mus=1/3*ones(nel,1)*0.3;
ref=1*ones(nel,1);
S_old=toastSysmat(hMesh,mua,mus,ref,0,'EL','PDD');

prm = ones(nel,1);
S=toastSysmatComponent(hMesh,prm,'PDD','EL');

subplot(224); spy(S)
img= calc_jacobian_bkgnd(imdl);
SS= calc_system_mat(img);
subplot(223); spy( SS.E)

gnd_node =1;
nnodes= size(qvec,1);
idx= (1:nnodes)'; idx(gnd_node)=[];
phi = zeros(nnodes,1);
phi(idx)=S(idx,idx)\qvec(idx);
meas= fwd_solve(img);
vv= meas.volt;
