#!/usr/bin/perl 
# beginnings of a file to do release assemply
# currently this is done by Makefile steps1 .. 7
use strict;
use warnings;

my $P = $ENV{'HOME'} . "/eidors-release/htdocs";
my $SZLIM = 30;
my %bigfiles;
open F, "find $P -size +$SZLIM -and -not -iname \*html |" or die $!;
while (<F>) { 
    s/\s*$//; $bigfiles{$_}++;
}
close F;
print %bigfiles;

open F, "find $P -iname \-s*html |" or die $!;

__DATA__
	find $P -name \*-s.html -exec \
	perl -i -pe'BEGIN{' -e'print STDERR "PROC: $$ARGV[0]\n";' \
             -e"%bigf = qw{$$bigfiles};" \
             -e'for (Pbigf) {s{%5C}{/}g; $$s=$$t= $$_;' \
             -e  '$$t=~ s{^.*(data_contrib|tutorial|news_pics|examples)(/.*)}{$(WEBLINK)$$1$$2};'\
             -e  '$$s=~ s{^.*/(.*)}{$$1};'\
             -e  '$$m{$$s}= $$t;'\
          -e'} } my $$q= q{['\''"]}; $$n= q{^['\''"]};' \
          -e'next unless m{(href|src) \s* = \s* ($$q)(.*/)?([^\2]*?)(\2) }x;' \
          -e'next unless $$m{$$4};' \
          -e's{(href|src) \s* = \s* ($$q)(.*/)?([^\2]*?)(\2) }' \
          -e '{$$1 = "$$m{$$4}"}x;' \
        \{} \;



