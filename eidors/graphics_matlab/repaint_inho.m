function repaint_inho(mat,mat_ref,vtx,simp, thresh, clr_def);
%function repaint_inho(mat,mat_ref,vtx,simp, thresh);
%
%Repaints the simulated inhomogeneity according to the reference
%distribution. (Increase -> Red, Decrease -> Blue) 
%
%mat     = The simulated (targeted) distribution.
%mat_ref = The known initial (homogeneous) distribution.
%        = Override default unless 'use_global'
%vtx     = The vertices matrix.
%simp    = The simplices matrix.
%thresh  = Threshold to show imaged region (or [] for default)
%clr_def = Colour definitions val.calc_colours.field etc

% (C) 2005 Andy Adler + Nick Polydorides. License: GPL version 2 or version 3
% $Id$

if nargin<5
    thresh = [];
end
if nargin<6
    clr_def = [];
end
if strcmp(mat_ref, 'use_global')
   img.calc_colours.ref_level = mat_ref;
end

if isempty(thresh)
    thresh = 1/4;
end

% looks best if eidors_colours.greylev < 0
[colours,scl_data] = calc_colours( mat, clr_def, 0);
ii=find( abs(scl_data) > thresh);
this_x = simp(ii,:);

colours= permute(colours(ii,:,:),[2,1,3]);
ELEM= vtx';

Xs=   zeros(3,length(ii));
Ys=   zeros(3,length(ii));
Zs=   zeros(3,length(ii));
switch(size(this_x,2))
    case 3
        idx_ = [1;2;3];
    case 4
        idx_ = [[1;2;3], ...
                [1;2;4], ...
                [1;3;4], ...
                [2;3;4]];
end
for idx=idx_
   Xs(:)=vtx(this_x(:,idx)',1);
   Ys(:)=vtx(this_x(:,idx)',2);
   Zs(:)=vtx(this_x(:,idx)',3);

   if exist('OCTAVE_VERSION');
% TODO: This is really slow, can we do anything about it
      cmap = colormap;
      for i=1:size(colours,2);
         patch(Xs(:,i),Ys(:,i),Zs(:,i),cmap(colours(i),:));
      end
   else
   if size(colours,1)==1 && size(colours,2)==3
      % need to work around ^%$#%$# matlab bug which
      % forces an incorrect interpretation is colours of this size
      hh= patch(Xs(:,[1:3,1]), ...
                Ys(:,[1:3,1]), ...
                Zs(:,[1:3,1]), ...
                colours(:,[1:3,1]), ...
            'EdgeColor','none');
   else
      hh= patch(Xs,Ys,Zs,colours, ...
            'EdgeColor','none');
   end
   end
end
