# meshgrid

this folder generates the meshgrid for you which you will eventually
use in raydata_raymatrix and finally in your inversion.

you have to change values at the top of

* makemymesh.m

accordingly. have a look at **marias_mesher.m** where I tried to clean
the code a bit up.

If you want a regular grit you don't need to provide a station list,
but if you wish to adjust your grid in areas where your station density
is increased, you have to create a station list in this shape:

**stla	stlo	net.stat.loc.chan	stel	stdp**

-20.4742              24.5139    XK.B04KH..BHZ           921.0        0.0

-0.1937             36.08294    1C.MCN1..HHZ          1872.0        0.0

be aware that the station list better should not have a header otherwise
matlab has some problems with it.

It's fine to stop the code, at the end of **makemymesh.m** are instructions
how to proceed.

## Needs fixing!!
Just in case: qhull probably will not work, just go to the terminal and
execute this one line there and continue in matlab.


