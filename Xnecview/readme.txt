X N E C V I E W
---------------
A program for visualizing NEC2 input and output data (Unix/X-windows).

Copyright (C) 1998-2011, Pieter-Tjerk de Boer

-----------------------------------------------------------------------
This program is free software; you can redistribute it and/or modify
it under the terms of version 2 of the GNU General Public License as
published by the Free Software Foundation.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
-----------------------------------------------------------------------

COMPILING:
----------
This should be trivial on Linux systems: just type 'make' in the
directory where the source has been unpacked.
I also supplied an Imakefile, which might make compilation on other
systems easier: typing 'xmkmf' will generate a Makefile, which can then
be used by 'make'.
Furthermore, you might want to change a few of the #define lines in config.h,
to suit your tastes; see the comments.

If the PNG library is not installed on your system, you should either
- use the Makefile after deleting the 'PNG=yes' line in it
or
- use xmkmf after doing   mv Imakefile_nopng Imakefile
Obviously, a thus-compiled version of xnecview cannot export to PNG files.

If you get an error about missing getopt.h while compiling, put a line
  #define NO_GETOPT
at the beginning of xnecview.c . You'll lose the support for command line
options, though.
This apparently happens if you have a non-GNU (compatible) version of getopt.

And if you get an arror about drem(), uncomment the last line in xnecview.h
and recompile.

Last but not least, note that xnecview needs version 2.0 or later of the
GTK+ libraries.

USAGE:
------
Just invoke xnecview with, as its command-line arguments, one or more names
of files containing NEC2 input (structure of an antenna model) and/or output
data (impedance, radiation pattern). Then xnecview will try to graphically
display whatever data is in those files, in one or two windows.
Window 1 shows a 3D plot of the structure (wires etc.) of the antenna,
and/or the spatial distribution of the radiation, and/or distribution of
the current and on the antenna elements and near electric and magnetic
fields.
Window 2 shows a set of graphs of several quantities (SWR, gain, etc.)
as a function of frequency.
Keyboard and mouse can be used to manipulate these views and export them
as postscript files.
For more details, see the manpage (i.e., the file xnecview.1x).

WHAT'S NEW:
-----------
This is version 1.36; changes compared to version 1.35 are only a few
minor bugfixes (see the file HISTORY) and support for NEC output files
that do not have impedance data (practically speaking: structures that
are excited by an incoming field).

HELP NEEDED:
------------
I'd like some help with the following issues:
- GH-cards (helix/spiral) with helix-length = 0. I don't understand these,
  and neither does my copy of NEC. According to doc on the WWW, this is a
  non-official extension, but I'd still like to support it.
- What should be the 0 dB gain reference in the gain plots if gain data
  is available at more than one frequency? At present, the maximum gain
  observed over all frequencies is used for this; this can be handy for
  comparing gain at different frequencies, but not necessarily for
  judging the radiation pattern at one particular frequency, due to the
  non-linearity of some of the gain scales.
- In a similar vain: should the currents distribution use the same scaling
  at all frequencies if currents are available at several frequencies?
  At present, the scaling is performed separately for every frequency.
- As should be clear from the man page, the interpretation of the currents
  display is somewhat non-trivial. Hints on how to make this clearer (either
  in the documentation, or by changes to the program) would be welcome.
- The hidden-line removal code is somewhat fragile; it needs to make some
  assumptions about the structure of the radiation pattern data. The code
  tries to check that the data satisfies all conditions, and reverts to
  wire-grid plotting if needed.
  I'd be very interested in any reports about observations of unexpected or
  incorrect behaviour; e.g., lines being visible that shouldn't be, or the
  other way around, or the code reverting to the wire mesh drawing where
  it would be desirable to still have hidden-line removal (so the code should
  be extended to handle that case).
- The number of variants of gain-vs-frequency plots is getting large, now
  that the polarization can also be selected; and still some useful
  quantities like the axial ratio are not available. Suggestions on
  making this more manageable are welcome. (One idea: rather than
  distinguishing between maxgain and vgain, there could be just one
  gain plot, and a separate selector to choose between max-gain,
  viewer-gain, and possibly others.)

Contributions, comments and bug reports are welcome, please send them to
pa3fwm@amsat.org .
