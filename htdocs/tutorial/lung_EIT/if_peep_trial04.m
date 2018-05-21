% Display Times courses $Id$

for loop = 1:2;
   if loop == 1; img = i_injury;
   else          img = i_treat;
   end

   subplot(2,2,loop);
   time = (0:size(img.elem_data,2)-1)/13; % Frame rate = 13/s

   raster= calc_slices( img );
   ROIs  = raster(ylocn, xlocn, :);
   ROIs  = permute(ROIs, [3,1,2]);

   %Normalize to its maximum
   for i=1:length(ylocn)
      ROIs(:,i) = - ROIs(:,i) / max(abs(ROIs(:,i)));
   end
   plot(time, ROIs);
   axis([0, max(time), -0.1, 1]);
   legend('1','2','3','4','Location','NorthWest')
   xlabel('time (s)')
   ylabel('normalized \Delta Z')
end

print_convert if_peep_trial04.png '-density 175'
