Data collected Feb 26, 2008.

Data was averaged since the images are static (no movement during measurements).

matlab script:
  % Load measured data
  for i=20:28;
    vv=eval(sprintf('eidors_readdata(''phantom_data/data/999999100%d.get'');',i));
    eval(sprintf('phantomdata_%d=mean(vv,2);',i));
  end
  save XYZ phantomdata_*;

The averaged data is stored in variables phantomdata_* where * is:
  % 20 no compression, no target
  % 21      ''       , glass target 
  % 22      ''       , metal target

  % 23 2 pt compression, no target
  % 24      ''       , glass target 
  % 25      ''       , metal target

  % 26 3 pt compression, no target
  % 27      ''       , glass target 
  % 28      ''       , metal target


