function data= find_artefacts( data)
% Copyright C. Gomez-Laberge, August 2009.
% $Id$
   for i=1:208
      dc = data(i,:);
      if norm(dc) < 1e-4; continue; end

      qq= quartiles( dc );
      ff= find( abs(dc - qq(2) )  > 4*(qq(3)-qq(1)) );
      if ff
        fprintf('Artefact Detected in Measurement %d\n', i);
        dc(ff) = qq(2);
        data(i,:) = dc;
     end 
   end