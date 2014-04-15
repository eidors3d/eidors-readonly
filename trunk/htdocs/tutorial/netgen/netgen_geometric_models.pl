my $FILE= '../../../eidors/meshing/netgen/ng_mk_geometric_models.m';
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
Cylinder with toruses (tori??)
</h3>

<pre>
LIST1

my $LIST2=<<'LIST2';
</pre>

<center>
<img src="netgen_cyl_models%02d.png">
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
