% Test algorithm performance $Id$

% Reconstruct GREIT Images
imdl_gr = mk_common_gridmdl('GREITc1');

% Reconstruct backprojection Images
imdl_bp = mk_common_gridmdl('backproj');

% Reconstruct GN Images
imdl_gn = select_imdl( mk_common_model('d2c2', 16), {'Basic GN dif','Choose NF=0.5'});

test_performance( { imdl_gr, imdl_bp, imdl_gn } );

print_convert 'algorithm_performance01a.png' '-density 100'
