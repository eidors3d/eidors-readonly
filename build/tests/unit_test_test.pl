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
print F "TEST= TEST + 1; fprintf('TEST=%d\\n', TEST);\n";
print F "switch TEST\n"; my $i=1;
foreach (sort keys %tuts) {
  print F "case $i; $_('UNIT_TEST')\n";
  $i++;
}
print F "end;\n";
close F;


