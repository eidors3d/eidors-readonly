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
print F "unit_test_cmp('RESET_COUNTER');\n";
print F "fprintf('TEST=%d fn=%s\\n', TEST, func2str(fn));\n";
print F "eidors_cache clear; clf; feval(fn, 'UNIT_TEST');\n";
print F "unit_test_cmp('SHOW_COUNTER',func2str(fn));\n";
close F;


