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
       m{^.*/(.*).m$}; $octcmd .= " $1;";
    }
    $octcmd .="'";
   print $octcmd;
   
   system $octcmd; 
last
}


__DATA__
my @files = File::Find::Rule->file()
              ->name( "*.m" ) ->in( "../../htdocs/tutorial" );

my %tuts;
for my $file (@files) {
   $tuts{$1}++ if $file =~ /(.*)[0-9][0-9].*\.m/;
}
  
foreach (keys %tuts) {
    my $files = `ls "$_*"`; 
    print "$_\n";
}
find ../../htdocs/tutorial/ -name *.m | sort | \
perl -ne'$tuts{$1}++ if /(.*)[0-9][0-9].*\.m/;' \
   -e'END{foreach (keys %tuts){' \
   -e' $files = `ls "$_*"`;' \
   -e' print "$_\n";'  \
   -e'};}'

 # `ls $_*.m`;'
