#!perl -w
use strict;

open F, "find ../../htdocs/tutorial/ -name *.m | sort |" or die $!;
my %tuts;
while ( <F> ) {
   $tuts{$1}++ if $_ =~ /(.*)[0-9][0-9].*\.m/;
}
close F;

foreach (sort keys %tuts) {
    my $cmd = "$_" . "*.m";
    my $octcmd = "octave -q --eval 'run /home/adler/docs/eidors/eidors/startup.m;";
    m{(.*)/[^/.*]}; $octcmd .= " cd $1; ";
    open F, "ls $cmd| sort|" or die $!;
    while ( <F> ) {
       print $_;
       m{^.*/(.*).m$}; $octcmd .= "disp([\"############## $1 ############\"]); $1;";
    }
    $octcmd .="'";
   print $octcmd;
   
   system "$octcmd";



