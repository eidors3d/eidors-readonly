[stim,mpat]= mk_stim_patterns(16,1,[0,1],[0,1],{},1);

[vh, vi, p] = face('large');

%vi = add_noise(4,vi,vh);
vi.meas(1:10:30) = 5;

%fmdl = ng_mk_cyl_models([0,1,.1],16,.05);
imdl = mk_common_model('c2c2',16);fmdl = imdl.fwd_model;

fmdl.stimulation = stim;
imdl = select_imdl(fmdl, {'Basic GN dif'});
imdl.solve = @inv_solve_diff_pdipm;
imdl = rmfield(imdl,'RtR_prior');
imdl.R_prior = @prior_TV;

for i=1:2; for j=1:2;
    switch 10*i+j;
        case 11; hp= .01;
        case 12; hp= .0001;
        case 21; hp= .1;
        case 22; hp= .01;
    end
    imdl.hyperparameter.value = 10*hp;
    
    imdl.inv_solve_diff_pdipm.norm_image = i;
    imdl.inv_solve_diff_pdipm.norm_data  = j;

    subplot(2,2,2*i + j - 2 )
    show_fem(inv_solve(imdl, vh,vi));
    title(sprintf('im = %d  data = %d',i,j));
end;end