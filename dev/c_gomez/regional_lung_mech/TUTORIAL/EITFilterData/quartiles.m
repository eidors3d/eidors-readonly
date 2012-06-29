function qq= quartiles( dd )
% Copyright C. Gomez-Laberge, August 2009.
% $Id$
   ds= sort( dd );
   qq= interp1( linspace(0,1,length(dd)), ds, [.25,.5,.75] ); 
 