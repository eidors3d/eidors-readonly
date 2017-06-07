#!perl -w
use strict;

open F, 'find ../../htdocs/tutorial/ -name \*.m | sort |' or die $!;
my @tuts;
while ( <F> ) {
   push( @tuts, $1) if $_ =~ m{(.*/.*\.m)};
}
close F;


open F, "> run_currenT.m" or die $!;
print F "if ~exist('TEST'); TEST =0; end; TEST= TEST + 1;\n";
print F "switch TEST\n"; my $i=1;
foreach (@tuts) {
  printf F "case %3d; fn='%s';\n", $i, $_;
  $i++;
}
print F "end;\n";
print F <<TESTCODE;
diary tutorial_test_out.txt
fprintf('TEST=%d fn=%s\\n', TEST, fn);
try
   run(fn);
catch
   L = lasterror();
   fprintf(' ERROR = (%s)\\n',L.message);
end
diary off
TESTCODE
close F;


