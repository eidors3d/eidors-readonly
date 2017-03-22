% create reconstruction models and select hyperparameter with different approaches
% $Id$
if ~exist('USE_GREIT_NOT_GN', 'var') || ~USE_GREIT_NOT_GN
    lambda_sel_approach = {'NF', 'SNR', 'LCC', 'GCV'};
    imdl_creation_fun = @(img, opts, lambda) mk_GN_model(img, opts, lambda);    
    USE_GREIT_NOT_GN = false;
else
    lambda_sel_approach = {'NF', 'SNR'};
    imdl_creation_fun = @(img, opts, lambda) mk_GREIT_model(img, 0.2, lambda, opts);
end
    
clf;
for jj = 1:length(lambda_sel_approach)
    for ii = 1:length(imgs)
        subplot(length(lambda_sel_approach), 4, ii + (jj-1)*length(imgs));
        
        opts = [];
        opts.img_size = [32 32];
        switch(lambda_sel_approach{jj})
            case 'NF'   % fixed noise figure
                lambda = [];    
                opts.noise_figure = 0.5;    
            case 'SNR'  % fixed image SNR as suggested by Braun et al. 
                lambda = 1;  % lambda is used as initial weight to find the appropriate SNR
                             % this is an arbitrary initial guess to avoid non-convergence
                opts.image_SNR = imdl{1,1}.SNR;   % same as NF=0.5 for 16 elecs adjacent
            case 'LCC'  % LCC: L-curve criterion
                lambda = 1; % lambda will be selected further below
            case 'GCV'  % GCV: generalized-cross validation
                lambda = 1; % lambda will be selected further below
            otherwise
                error(['Unknown hyperparameter selection approach: ', lambda_sel_approach(jj)]);
        end
        
        % create reconstruction framework
        imdl{jj, ii} = imdl_creation_fun(mk_image(rmdls{ii},1), opts, lambda);
        imdl{jj, ii}.fwd_model.name = imgs{ii}.fwd_model.name;
        
        % add noise to inhomogeneous conductivity change
        vn = repmat(vi{ii}.meas, 1, size(noise,2)) + noise(1:length(vi{ii}.meas), :);
        
        if strcmp(lambda_sel_approach{jj}, 'LCC') || strcmp(lambda_sel_approach{jj}, 'GCV')
            % use simulated noisy data to choose hyperparameter either with:
            % LCC: L-curve criterion
            % GCV: generalized-cross validation
            lambdas = calc_lambda_regtools(imdl{jj, ii}, vh{ii}.meas, vn, lambda_sel_approach{jj});
            imdl{jj, ii}.hyperparameter.value = median(lambdas);
        end
        
        % calulate image SNR for each approach
        imdl{jj, ii}.SNR = calc_image_SNR(imdl{jj, ii});
        
        % perform image reconstruction
        imgr{jj, ii} = inv_solve(imdl{jj, ii}, vh{ii}.meas, vn);
                
        % calulate temporal RMS image as in Braun et al., 2017, IEEE TBME        
        trmsa{jj, ii} = imgr{jj, ii};
        trmsa{jj, ii}.elem_data = rms(trmsa{jj, ii}.elem_data, 2);
        % normalize to maximal amplitude
        norm_value = max(trmsa{jj, ii}.elem_data);
        trmsa{jj, ii}.elem_data = trmsa{jj, ii}.elem_data / norm_value;
        imgr{jj, ii}.elem_data = imgr{jj, ii}.elem_data / norm_value;
        % force displaying values in range [0 1]
        trmsa{jj, ii}.calc_colours.clim = 0.5;
        trmsa{jj, ii}.calc_colours.ref_level = 0.5;
        trmsa{jj, ii}.calc_colours.cmap_type = 'greyscale';
        imgr{jj, ii}.calc_colours.clim = 1;
        imgr{jj, ii}.calc_colours.ref_level = 0;
        
        % show tRMSA image for each configuration
        h = show_fem(trmsa{jj, ii}, 1);
        set(h, 'edgecolor', 'none');        
        title([num2str(ii, '(%01d)'), ' ', lambda_sel_approach{jj}, ': ', imdl{jj, ii}.fwd_model.name]);
        set(gca, 'ytick', [], 'xtick', []);
        if ii == 1
            ylabel({['(', char(double('a')+jj-1), ') ', lambda_sel_approach{jj}], 'Right'});
        else
            ylabel('Right');
        end
        xlabel({'Dorsal', sprintf('SNR = %03.2d, \\lambda = %03.2d', imdl{jj, ii}.SNR, imdl{jj, ii}.hyperparameter.value)});
    end
end
    
set(gcf, 'paperpositionmode', 'auto');
if USE_GREIT_NOT_GN
    print_convert np_trmsa_GREIT.png '-density 100'
else
    print_convert np_trmsa_GN.png '-density 100'    
end
