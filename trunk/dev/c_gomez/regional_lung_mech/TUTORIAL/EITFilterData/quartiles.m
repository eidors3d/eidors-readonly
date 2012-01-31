function qq= quartiles( dd )
   ds= sort( dd );
   qq= interp1( linspace(0,1,length(dd)), ds, [.25,.5,.75] ); 
 