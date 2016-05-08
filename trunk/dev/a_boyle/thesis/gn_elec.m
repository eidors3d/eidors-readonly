% (C) 2016, Alistair Boyle
% License: GPL version 2 or version 3

Wm=W(vb_recip,0.001);%*spdiag(Wn1);
f_mv = imdlA.fwd_model.mv_elec;
a =  [ 28.3393   43.5668   68.0494   79.9218 ]; % zoom in on elec#9 movement
eidors_msg log_level 2;
clear img vv str;
img = imgAa;
img.theta = theta;
img.elem_data = 1./imgAa.elem_data;
img = rmfield(img,'current_params');
img(1) = img;
imgH = img(1); imgH.elem_data = ones(size(imgH.elem_data)); vh = fwd_solve(imgH); vh = vh.meas;
if 1; R = eye(32); % Tikhonov
else; R = -spdiag(ones(33,1),-1) + spdiag(ones(32,1)*2,0) - spdiag(ones(31,1),1); end % Laplace 1D smoothing
if 1;             N = @(dx) dx;            Ninv = @(dx) dx;             dN = @(dx) 1; % direct movement
else; delta=1e-0; N = @(dx) log(dx+delta); Ninv = @(dx) (10.^dx)-delta; dN = @(dx) spdiag(N(sum(dx,2)+delta)); end % log movment
if 1; G = @(dv,vh) dv;     dG = @(dv,vh) 1;                 % measured voltages
else; G = @(dv,vh) dv./vh; dG = @(dv,vh) spdiag(1./vh); end % apparent resistivity
%if 1; fixed_elec = [1:3,20,30:32]; % don't allow electrodes at the ends to move
if 1; fixed_elec = [1:3,30:32]; % don't allow electrodes at the ends to move
%if 1; fixed_elec = [1:6]; % don't allow electrodes at the downslope end to move
else; fixed_elec = []; end
if sel_line == 1
   hp = [1e-1 2e-2]; thres = 1e-2; max_dx = 1.6; % config for line#1
else
   %hp = [1e-1 5e-2 4e-2 3e-2 2e-2]; hp=hp(:); thres = 1e-4; max_dx = 1.60; % config for line#5
   hp = [ 7e-2 5e-2 3e-2 1e-2]; hp=hp(:); thres = 1e-2; max_dx = 1.60; % config for line#5
end
dx = zeros(32,1); r = 1e300*[1 1]; hp2RtR = hp(1)^2*R'*R; hpsel = 1; % init
dx(fixed_elec) = duvw(fixed_elec,2)-0.04; % CHEAT: give the fixed electrodes the correct locations!
x0 = dx;
Jmn = NaN; % last Jacobian
for i = 1:20
   vvt = fwd_solve(img(i)); vv(:,i) = vvt.meas;
   if 1; vat = va; vv1 = vv(:,1); % difference
   else; vat = 0; vv1 = 0; end    % absolute
   dv = G((vb-vat) - (vv(:,i)-vv1), vh);
   x = sum(dx,2);
   de = x - x0;
   r(i+1,:) = [dv'*Wm*dv   de'*hp2RtR*de]/2;
   drp = (sum(r(i+1,:) - r(i,:)))/sum(r(i,:));
   fprintf('   r=%0.9g = %0.9g + %0.9g (%0.7g%%)\n', sum(r(i+1,:)), r(i+1,1), r(i+1,2), drp*100);
   if drp > -thres; if hpsel == length(hp); disp('STOP'); break; else hpsel=hpsel+1; fprintf('   ** GEAR SHIFT ** hp=%g->%g\n',hp(hpsel-1),hp(hpsel)); end; end
   imgJ = img(i);
   if 1 imgJ.elem_data = mean(imgJ.elem_data)*ones(size(imgJ.elem_data)); end % do movement Jacobian on homogeneous estimate
   Jm = dG(dv,vh) * jacobian_movement_only(imgJ) * dN(dx);
   Jm = - Jm; % TODO Jm is negated... something stupid is happening
   Jm(:,fixed_elec) = 0; % force zero movement
   hpi = hp(hpsel); hp2RtR = hpi^2*R'*R;
   dx(:,i+1) = -(Jm'*Wm*Jm + hp2RtR)\(Jm'*Wm*dv - hp2RtR*(sum(dx,2)-x0));
   fprintf('iter#%d: hp=%0.2g, ||dJ||=%0.2g, ||dx||=%0.2g, max(dx)=%0.2g, min(dx)=%0.2g\n',i, hpi, norm(full(Jm-Jmn)), norm(dx(:,i+1)), max(dx(:,i+1)), min(dx(:,i+1))); Jmn = Jm;
%  if 1; dx(:,i+1) = linesearch(img(1), dx(:,i+1), vb, vat, vv1, hp2RtR, Jm, Wm, x, x0, G, Ninv, vh, f_mv, max_dx/max(abs(dx(:,i+1))), thres) * dx(:,i+1); end % line search // TODO fails with N != 1
   if 1; dx(:,i+1) = linesearch(img(1), dx(:,i+1), vb, vat, vv1, hp2RtR, Jm, Wm, x, x0, G, Ninv, vh, f_mv, 1, thres) * dx(:,i+1); end % line search // TODO fails with N != 1
   img(i+1) = img(i); img(i+1).fwd_model = f_mv(img(1).fwd_model, (Ninv(sum(dx,2))*[1 1]).*[cos(theta') sin(theta')] );
   for j=1:i+1; str{j} = sprintf('iter#%d',j-1); end
   ddu = [NaN; diff(duvw(:,2))]; sddx = [NaN; diff(sum(dx,2))];
%   clf; s0=svd(full(Jm'*Wm*Jm)); s1=svd(full(Jm'*Wm*Jm+hpi^2*R'*R)); plot([s0 s1]); legend('J^TWJ','J^TWJ+\lambda^2R^TR');title(sprintf('\\lambda=%g',hpi)); drawnow; pause(1);
   clf; bar([ duvw(:,2)-0.04 ddu sddx Ninv(cumsum(dx,2)) abs(duvw(:,2)-0.04-sum(dx,2)) abs(ddu-sddx)]); hold on; plot([1 32],[0.2 0.2],'--'); hold off; legend({ 'expected', 'expected spacing', 'spacing', str{:}, '|error|','|spacing error|'},'Location','EastOutside'); drawnow;
%   clf; show_fem(img(i+1),[0 1]); axis(a); drawnow;
end

   % we can cross check the internal calculations of linesearch() by looking at
   % the results returned by the anonymous functions ff{1:5}, they should
   % ultimately agree with our outer residual calculator
   idx = min(i+1,size(dx,2));
   [~,ff] = linesearch(img(1), dx(:,idx), vb, vat, vv1, hp2RtR, Jm, Wm, x, x0, G, Ninv, vh, f_mv, max_dx/max(abs(dx(:,idx))), 0, -1); % line search functions

   % THESIS figure...
common_ylim = [-1.8 +0.7];
   clf; show_fem(imgAa,1); % fix colors that got mucked up?
   clf; bar([  Ninv(cumsum(dx,2)) ]); legend({ str{:}},'Location','EastOutside');
        axis tight; box off; xlabel('(electrode #, iteration #)'); ylabel('electrode movement [m]'); drawnow;
   print('-dpdf', ['results/' name '-imgAMbi.pdf']);
   eval(['! pdfcrop results/' name '-imgAMbi.pdf']);
   eval(['! mv results/' name '-imgAMbi-crop.pdf results/' name '-imgAMbi.pdf']);

   % THESIS figure!
   clf; subplot(211); h=bar([ duvw(:,2)-0.04 (sum(dx,2))], 'EdgeColor', 'none');
        %hold on; plot([1:32]'*[1 1],[0.2*ones(1,32)]'*[1 -1] + duvw(:,2)*[1 1],'--','Color',[0.5 0.5 0.5]); hold off;
        hold on; plot([1:32]'*[1 1],[0.2*ones(1,32)]'*[1 -1],'--','Color',[0.5 0.5 0.5]); hold off;
        legend({'true movement','reconstructed'},'Location','SouthEast','Orientation','Horizontal','Box','off'); box off;
        ylabel('electrode mvmt [m]'); xlabel('electrode #'); axis tight;
        set(h(1),'FaceColor',[0.85 0.33 0.10]);
        set(h(2),'FaceColor',[0.47 0.67 0.19]);
        text(20,+0.35,'upslope'); text(20,-0.35,'downslope'); text(25,0.35,'0.2 m');
        ylim(common_ylim);
   print('-dpdf', ['results/' name '-imgAMb.pdf']);
   eval(['! pdfcrop results/' name '-imgAMb.pdf']);
   eval(['! mv results/' name '-imgAMb-crop.pdf results/' name '-imgAMb.pdf']);

   % THESIS figure!
   clf; subplot(211); h=bar(duvw,'EdgeColor','none'); box off;
        ylabel('electrode mvmt [m]'); xlabel('electrode #'); axis tight;
        set(h(1),'FaceColor','flat');
        set(h(2),'FaceColor',[0.85 0.33 0.10]);
        set(h(3),'FaceColor',[0.00 0.57 0.70]);
        ll=legend({'transverse','logitudinal','normal'},'Position',[0.6 0.65 0 0],'Orientation','Horizontal','Box','off');
        ylim(common_ylim);
   print('-dpdf', ['results/' name '-imgt.pdf']);
   eval(['! pdfcrop results/' name '-imgt.pdf']);
   eval(['! mv results/' name '-imgt-crop.pdf results/' name '-imgt.pdf']);
