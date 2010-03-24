% Simple model with three electrodes $Id$
elec_pts= {};
n_elecs= 8;
elec_width= 0.1;

hw= elec_width/2;
th = linspace(0,2*pi,n_elecs+1); th(end)=[];
for i=1:n_elecs;
   ti = th(i) + [hw;-hw];
   elec_pts{i} = [sin(ti),cos(ti)];
end

params = [0.10,10,0.05;
          0.10,10,0.02;
          0.10,10,0.005;
          0.10,30,0.005;
          0.10,50,0.005;
          0.05,10,0.02;
          0.03,10,0.02;
          0.02,10,0.02;
          0.10,10,0.001;
          0.05,10,0.001 ];

clf
for i=1:5; for j=0:1
   param= params(i+5*j,:);
   fmdl= dm_2d_circ_pt_elecs( elec_pts, [], param );

   axes('position',[0.05+0.49*j,(i-1)*0.17+0.10,0.45,0.13]);
   show_fem(fmdl);
   axis([0,1.05,-0.2,0.2]);
   text(0,0.25,sprintf('Params=[%4.2f,%2d,%5.3f]. #Elems= %d ', param, size(fmdl.elems,1)));
   if i>1; set(gca,'XTickLabel',[]); end
end; end

print -dpng -r200 circ_2d_model03.png
