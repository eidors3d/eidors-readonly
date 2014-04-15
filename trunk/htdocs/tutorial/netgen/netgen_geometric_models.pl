my ($FOUT1, $FOUT2, $FIN, $HEAD, $LIST1, $LIST2, $TAIL, $MCODE) = names();
open FOUT1,">$FOUT1" or die "Can't open $FOUT1: $!";
open FOUT2,">$FOUT2" or die "Can't open $FOUT2: $!";
open FIN,  "<$FIN"   or die "Can't open $FIN  : $!";
print FOUT2 $HEAD;

my $seen_do_test_number = 0;
my $start_case = 0;
my $this_case = -1; # not seen yet
my $new_case  =  0;
my $description = '';
my $case_text = '';
my $cut_blanks = '';
while( <FIN> ) {
   next if $this_case ==0; # last marker
   $seen_do_test_number = 1 if /function fmdl = do_test_number/;
   next unless $seen_do_test_number;
   #printf "%2d, %2d, %2d, %2d:  %s", $this_case, $start_case, length($description),length($cut_blanks), $_;
   if( /^ *case (\d+)/ ) {
     $new_case = $1;
     $start_case = 0;
     if ($case_text) {
        printf FOUT2 $LIST1, $this_case, $description;
        print  FOUT2 $case_text;
        printf FOUT2 $LIST2, $this_case;

        printf FOUT1 "\n\n%CASE %2d %%%%%%%%\n", $this_case;
        printf FOUT1 "disp('#### %02d ####');clear;\n", $this_case;
        print  FOUT1 $case_text;
        printf FOUT1 $MCODE, $this_case;
        $case_text = '';
     }
     $this_case = $new_case;
     next;
   }
   if( /^( *)%DESC: (.*)/ ) {
     $start_case = 1;
     $cut_blanks = $1;
     $description = $2;
     $description =~ s/\r$//;
     next;
   }
   next if $start_case == 0;
   if ($start_case) {
      my $line = $_;
      $line =~ s/^$cut_blanks//;
      $line =~ s/\r$//;
      $case_text .= $line;
   }
}
close FIN;

print FOUT2 $TAIL;
close FOUT2;


sub names  {
my $FOUT1= 'netgen_geometric_models.m';
my $FOUT2= 'netgen_geometric_models.shtml';
my $FIN= '../../../eidors/meshing/netgen/ng_mk_geometric_models.m';
my $HEAD=<<'HEAD';
<!--#set var="root" value="../../" -->
<!--#set var="show_tut" value="1" -->
<!--#include virtual="../../nav-sidebar.shtml" -->

<h2> 
Using new <tt>ng_mk_geometric_models</tt> interface
</h2>

EIDORS can use
<a href="http://sourceforge.net/projects/netgen-mesher/">
Netgen</a> to create sophisticated 2D and 3D models

<p>
Here are some examples of the varity of models which
can be generated using the function: <tt>ng_mk_geometric_models</tt>.
<p>
HEAD

my $LIST1=<<'LIST1';
<h3>
%d: %s
</h3>

<pre>
LIST1

my $LIST2=<<'LIST2';
</pre>

<center>
<img src="netgen_geometric_models%02d.png">
</center>
LIST2

my $TAIL=<<'TAIL';
</td></tr></table>
<p>
<small>
    Last Modified: $Date$ by $Author$
</small>
</BODY></HTML>
TAIL

my $MCODE=<<'MCODE';
show_fem( fmdl );
print_convert netgen_geometric_models%02d.png
MCODE

return ($FOUT1, $FOUT2, $FIN, $HEAD, $LIST1, $LIST2, $TAIL, $MCODE);
}
