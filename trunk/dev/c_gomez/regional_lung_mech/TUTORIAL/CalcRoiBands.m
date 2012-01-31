function [bands] = CalcRoiBands(lungROI)
%CALCROIBANDS   Describe function here.
%   
% Signature:
% function [out1,out2] = FunctionTemplate(in1,in2,in3)
%
% Input:
% in1       double      scalar      description here
% in2       double      MxN         description here
% (in3)     boolean     scalar      description here
% 
% Output:   
% out1      double      scalar      description here
% out2      boolean     scalar      description here
% s         boolean     scalar      errors present: true
% 
% Copyright C. Gomez-Laberge, June 2011.
% $Id: $

vdproj = sum(lungROI,2);

m50 = sum(vdproj.*[1:32]')/sum(vdproj);
floorm50 = floor(m50);
m25 = sum(vdproj(1:floorm50).*[1:floorm50]')/sum(vdproj(1:floorm50));
ceilm50 = ceil(m50);
m75 = sum(vdproj(ceilm50:32).*[ceilm50:32]')/sum(vdproj(ceilm50:32));

canvas = zeros(size(lungROI));
ventdist = canvas;
ventdist(1:floor(m25),:)=1;
ventmed = canvas;
ventmed(ceil(m25):floor(m50),:)=1;
dorsmed = canvas;
dorsmed(ceil(m50):floor(m75),:)=1;
dorsdist = canvas;
dorsdist(ceil(m75):32,:)=1;
bands{1} = logical(ventdist.*lungROI);
bands{2} = logical(ventmed.*lungROI);
bands{3} = logical(dorsmed.*lungROI);
bands{4} = logical(dorsdist.*lungROI);

% End of function
end %function