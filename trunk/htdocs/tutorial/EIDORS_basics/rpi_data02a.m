% load RPI tank data $Id$

load Rensselaer_EIT_Phantom;
vh = real(ACT2006_homog);
vi = real(ACT2000_phant);

% Preprocessing data. We ensure that each voltage sequence sums to zero
  for i=0:0 %30
    idx = 32*i + (1:32);
    vh(idx) = vh(idx) - mean(vh(idx));
    vi(idx) = vi(idx) - mean(vi(idx));
  end

 load ~/docs/carleton/2010/projects/lionheart-EIT-book/pre_proc_data.mat

