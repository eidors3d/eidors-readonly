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

mua=zeros(nel,1);
qvec = toastQvec(hMesh,'Neumann','Point');
mus=ones(nel,1);
ref=ones(nel,1);
S=toastSysmat(hMesh,mua,mus,ref,0,'EL');
subplot(224); spy(S)
img= calc_jacobian_bkgnd(imdl);
SS= calc_system_mat(img);
subplot(223); spy( SS.E)

phi=S\qvec;
meas= fwd_solve(img);
vv= meas.volt;