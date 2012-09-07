function stim= stim_pattern_geophys( n_elec, pat_type,  options )
% STIM_PATTERN_GEOPHYS: Create Geophysical Stimulation Patterns
%
% stim= stim_pattern_geophys( n_elec, type, options_list )
%  n_elec = Number of electrodes (in single line)
%  type   = name of stimulation pattern (string)
%  options_list = options for the pattern
%
%  options = {'current', 0.1}  % 0.1 Amps current (default 1)
%  options = {'gain', 2.0}     % 2.0 gain         (default 1)
%  options = {'spacings',[1,2,3,4]} list of spacings 
%  options = {'reciprocal_meas',1} % do reciprocal (default 0);
%
% TYPE = 'WENNER'
% 

% (C) 2011 Andy Adler. License: GPL version 2 or version 3
% $Id$

if isstr(n_elec) && strcmp(n_elec,'UNIT_TEST'); do_unit_test; return; end
if nargin<3; options = {}; end
options = parse_options(options);

switch(upper(pat_type))
    case 'WENNER';          stim = stim_pattern_wenner(n_elec, options);
    case 'SCHLUMBERGER';    stim = stim_pattern_schlumberger(n_elec, options);
    case 'DIPOLEDIPOLE';   stim = stim_pattern_dipoledipole(n_elec, options);
  otherwise;
    error('No pattern of type "%s" available', upper(pat_type));
end

stim = stim_meas_list( stim, n_elec, options.current, options.gain); 

function pp= parse_options(opts)
% Set defaults here
%   pp.spacings = [1,2,4,8,16,32,64,128,256,1024,2048];
  pp.spacings= 1:256;
  pp.multiples= ones(1,256);
  pp.reciprocal_meas = 0;
  pp.current = 1;
  pp.gain = 1;

  for i= 1:2:length(opts)-1;
    pp.( lower( opts{i} ) ) = opts{i+1};
  end

function stim = stim_pattern_wenner(n_elec, options)
  stim=[];
  for a= options.spacings;
    block =  [0,3,1,2]*a + 1;
    if max(block)>n_elec
       break;
    end
    for i=(0:n_elec-max(block))
       stim = [stim; i+block];
       if options.reciprocal_meas
          stim = [stim; i+block([3,4,1,2])];
       end
    end
  end

function stim = stim_pattern_schlumberger(n_elec, options)
  stim=[];
  a= options.spacings;
  n= options.multiples;
  if length(a)~=length(n)
   error('spacings and multiples must be vector with a same size');
  end
  for k= 1:length(a);
      block =  [0,round((2*n(k)+1)*a(k)),round(n(k)*a(k)),round((n(k)+1)*a(k))]+1;
    if max(block)>n_elec
       break;
    end
    for i=(0:n_elec-max(block))
       stim = [stim; i+block];
       if options.reciprocal_meas
          stim = [stim; i+block([3,4,1,2])];
       end
    end
  end  
  
  
function stim = stim_pattern_dipoledipole(n_elec, options)
  stim=[];
  a= options.spacings;
  n= options.multiples;
  if length(a)~=length(n)
   error('spacings and multiples must be vector with a same size');
  end
  for k= 1:length(a);
      block =  [0,round((n(k)+2)*a(k)),a(k),round((n(k)+1)*a(k))]+1;
    if max(block)>n_elec
       break;
    end
    for i=(0:n_elec-max(block))
       stim = [stim; i+block];
       if options.reciprocal_meas
          stim = [stim; i+block([3,4,1,2])];
       end
    end
  end  
  
function do_unit_test
  unit_test_wenner;
  unit_test_schlumberger;
  unit_test_dipoledipole;

function unit_test_wenner
  stim= stim_pattern_wenner( 13, parse_options({}) ); 
  test1= [1,4,2,3; 2,5,3,4;
3,6,4,5; 4,7,5,6; 5,8,6,7; 6,9,7,8; 7,10,8,9; 8,11,9,10;
9,12,10,11; 10,13,11,12; 1,7,3,5; 2,8,4,6; 3,9,5,7;
4,10,6,8; 5,11,7,9; 6,12,8,10; 7,13,9,11; 1,10,4,7; 2 11 5 8; 3 12 6 9;
4 13 7 10; 1,13,5,9];
  unit_test_cmp('WENNER #1', stim, test1);
  stim= stim_pattern_geophys( 13,'wenner', {});
  unit_test_cmp('WENNER #2', stim, stim_meas_list(test1));


  stim= stim_pattern_geophys( 5,'wenner', {'reciprocal_meas',0});
  unit_test_cmp('WENNER #3a',stim, stim_meas_list([1,4,2,3;2,5,3,4]));
  stim= stim_pattern_geophys( 5,'wenner', {'reciprocal_meas',1});
  unit_test_cmp('WENNER #3b',stim, stim_meas_list([1,4,2,3;2,3,1,4;2,5,3,4;3,4,2,5]));

  stim= stim_pattern_geophys( 5,'wenner', {'current',0.1});
  unit_test_cmp('WENNER #4a',stim, stim_meas_list([1,4,2,3;2,5,3,4],5,0.1));
  stim= stim_pattern_geophys( 5,'wenner', {'gain',0.1});
  unit_test_cmp('WENNER #4b',stim, stim_meas_list([1,4,2,3;2,5,3,4],5,1,0.1));

  stim= stim_pattern_geophys( 5,'wenner', {'spacings',[1,2,3]});
  unit_test_cmp('WENNER #5a',stim, stim_meas_list([1,4,2,3;2,5,3,4]));
  stim= stim_pattern_geophys( 7,'wenner', {'spacings',[2,1,3]});
  unit_test_cmp('WENNER #5b',stim, stim_meas_list(...
       [1,7,3,5; 1,4,2,3; 2,5,3,4; 3,6,4,5; 4,7,5,6]));


function unit_test_schlumberger
    stim= stim_pattern_schlumberger( 13, parse_options({}) );
    test1= [1,4,2,3; 2,5,3,4;
        3,6,4,5; 4,7,5,6; 5,8,6,7; 6,9,7,8; 7,10,8,9; 8,11,9,10;
        9,12,10,11; 10,13,11,12; 1,7,3,5; 2,8,4,6; 3,9,5,7;
        4,10,6,8; 5,11,7,9; 6,12,8,10; 7,13,9,11; 1,10,4,7; 2 11 5 8; 3 12 6 9;
        4 13 7 10; 1,13,5,9];
    unit_test_cmp('SCHLUMBERGER #1', stim, test1);
    stim= stim_pattern_geophys( 13,'schlumberger', {});
    unit_test_cmp('SCHLUMBERGER #2', stim, stim_meas_list(test1));


    stim= stim_pattern_geophys( 5,'schlumberger', {'reciprocal_meas',0});
    unit_test_cmp('SCHLUMBERGER #3a',stim, stim_meas_list([1,4,2,3;2,5,3,4]));
    stim= stim_pattern_geophys( 5,'schlumberger', {'reciprocal_meas',1});
    unit_test_cmp('SCHLUMBERGER #3b',stim, stim_meas_list([1,4,2,3;2,3,1,4;2,5,3,4;3,4,2,5]));

    stim= stim_pattern_geophys( 5,'schlumberger', {'current',0.1});
    unit_test_cmp('SCHLUMBERGER #4a',stim, stim_meas_list([1,4,2,3;2,5,3,4],5,0.1));
    stim= stim_pattern_geophys( 5,'schlumberger', {'gain',0.1});
    unit_test_cmp('SCHLUMBERGER #4b',stim, stim_meas_list([1,4,2,3;2,5,3,4],5,1,0.1));

    stim= stim_pattern_geophys( 6,'schlumberger', {'spacings',[1,1],'multiples',[1,2]});
    unit_test_cmp('SCHLUMBERGER #5a',stim, stim_meas_list([1,4,2,3;2,5,3,4;3,6,4,5;1,6,3,4]));
    stim= stim_pattern_geophys( 11,'schlumberger', {'spacings',[1,2],'multiples',[3,2]});
    unit_test_cmp('SCHLUMBERGER #5b',stim, stim_meas_list(...
        [1,8,4,5;2,9,5,6;3,10,6,7;4,11,7,8;1,11,5,7]));
   
% options.spacing= [1 1 1 2 3 4 6 8 8 11 12 14 17];
% options.multiples= [1 2 3 2 5/3 6/4 7/6 1 10/8 1 13/12 15/14 1];

function unit_test_dipoledipole
    stim= stim_pattern_dipoledipole( 13, parse_options({}) );
    test1= [1,4,2,3; 2,5,3,4;
        3,6,4,5; 4,7,5,6; 5,8,6,7; 6,9,7,8; 7,10,8,9; 8,11,9,10;
        9,12,10,11; 10,13,11,12; 1,7,3,5; 2,8,4,6; 3,9,5,7;
        4,10,6,8; 5,11,7,9; 6,12,8,10; 7,13,9,11; 1,10,4,7; 2 11 5 8; 3 12 6 9;
        4 13 7 10; 1,13,5,9];
    unit_test_cmp('DIPOLEDIPOLE #1', stim, test1);
    stim= stim_pattern_geophys( 13,'dipoledipole', {});
    unit_test_cmp('DIPOLEDIPOLE #2', stim, stim_meas_list(test1));


    stim= stim_pattern_geophys( 5,'dipoledipole', {'reciprocal_meas',0});
    unit_test_cmp('DIPOLEDIPOLE #3a',stim, stim_meas_list([1,4,2,3;2,5,3,4]));
    stim= stim_pattern_geophys( 5,'dipoledipole', {'reciprocal_meas',1});
    unit_test_cmp('DIPOLEDIPOLE #3b',stim, stim_meas_list([1,4,2,3;2,3,1,4;2,5,3,4;3,4,2,5]));

    stim= stim_pattern_geophys( 5,'dipoledipole', {'current',0.1});
    unit_test_cmp('DIPOLEDIPOLE #4a',stim, stim_meas_list([1,4,2,3;2,5,3,4],5,0.1));
    stim= stim_pattern_geophys( 5,'dipoledipole', {'gain',0.1});
    unit_test_cmp('DIPOLEDIPOLE #4b',stim, stim_meas_list([1,4,2,3;2,5,3,4],5,1,0.1));

    stim= stim_pattern_geophys( 6,'dipoledipole', {'spacings',[1,1],'multiples',[1,2]});
    unit_test_cmp('DIPOLEDIPOLE #5a',stim, stim_meas_list([1,4,2,3;2,5,3,4;3,6,4,5;1,5,2,4;2,6,3,5]));
    stim= stim_pattern_geophys( 9,'dipoledipole', {'spacings',[1,2],'multiples',[3,2]});
    unit_test_cmp('DIPOLEDIPOLE #5b',stim, stim_meas_list(...
        [1,6,2,5;2,7,3,6;3,8,4,7;4,9,5,8;1,9,3,7]));
    
    
    