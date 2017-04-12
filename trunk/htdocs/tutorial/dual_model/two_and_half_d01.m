% Models
fmdl = mk_common_model('d2C',16); % fine model
cmdl = mk_common_model('c2C',16); % coarse model
c2f = mk_coarse_fine_mapping(fmdl.fwd_model, cmdl.fwd_model);
imdl = fmdl;
imdl.rec_model = cmdl.fwd_model;
imdl.fwd_model.coarse2fine = c2f;

clf; % figure 1
subplot(121);show_fem(imdl.fwd_model);
title('fine (2d) model'); axis square; axis off;
subplot(122); show_fem(imdl.rec_model);
title('coarse (2d) model'); axis square; axis off;
print_convert two_and_half_d01a.png '-density 75'
