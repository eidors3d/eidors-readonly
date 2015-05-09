#!perl -w
use strict;

open F, 'find ../../eidors -name \*.m | xargs grep -l UNIT_TEST | sort |' or die $!;
my %tuts;
while ( <F> ) {
   $tuts{$1}++ if $_ =~ m{.*/(.*)\.m};
}
close F;


open F, "> run_current.m" or die $!;
print F "if ~exist('TEST'); TEST =0; end;\n";
print F "TEST= TEST + 1;\n";
print F "switch TEST\n"; my $i=1;
foreach (sort keys %tuts) {
  print F "case $i; fn=\@$_;\n";
  $i++;
}
print F "end;\n";
print F <<TESTCODE;
diary unit_test_out.txt
unit_test_cmp('RESET_COUNTER');
fprintf('TEST=%d fn=%s\\n', TEST, func2str(fn));
eidors_cache clear; clf;  close all;
try
   feval(fn, 'UNIT_TEST');
   unit_test_cmp('SHOW_COUNTER',func2str(fn));
catch
   L = lasterror();
   fprintf(' ERROR = (%s)\\n',L.message);
end
diary off
TESTCODE
close F;


