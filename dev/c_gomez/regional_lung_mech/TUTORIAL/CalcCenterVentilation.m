function [cov,area] = CalcCenterVentilation(lungROI)
% For JK's manuscript on the RERV study November 2011

vdproj = sum(lungROI,2);
m50 = sum(vdproj.*[1:32]')/sum(vdproj);
cov = round(100*(m50-16)/16)/100;
area = round(100*sum(lungROI(:))/32^2);
