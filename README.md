Rivet Jet Substructure Study
============================

Authors
-------
Ayana Arce <atarce@phy.duke.edu>

David Bjergaard <david.b@duke.edu>

Deepak Kar <deepak.kar@cern.ch>

Karl Nordstrom <karl.nordstrom@cern.ch>

Rivet implementation
--------------------

The code is implemented as an additional rivet
projection called "BOOSTFastJets", to see how
to use it please have a look at the header
(in ./include) and the MC_GENSTUDY_JET*.cc analysis files.

Obtaining and running
---------------------

1. Setup the rivet installation you want to use with
(in the correct directory):

`source rivetenv.sh`

2. Patch FastJets to include the particles() method,
recompile FastJets:

`MY_RIVET=/path/to/your/rivet`

`THIS_DIR=$(pwd)`

`patch  $MY_RIVET/include/Rivet/Projections/FastJets.hh < $THIS_DIR/BOOSTFastJets.patch`

`cd $MY_RIVET && make -j 4 && make install`

3. While in this directory, run:

`make && make install` **

4. Have fun looking at substructure histograms!

** If you need super user privileges to do a make install
you'll have to sudo this; if your super user environment
doesn't have the correct rivet paths set up
(rivet-config: command not found) you'll have to do things
manually:

`LIBDIR=$(rivet-config --libdir)`

`sudo cp -f libBOOSTFastJets.so $LIBDIR`

Building an analysis and linking our library
--------------------------------------------

The makefile should provide an example of how to build an analysis
that uses code in the BOOSTFastJets projection. If you want to do this
yourself, keep in mind that you will have to explicitly link to the
library when using rivet-buildplugin (in other words, if you decide
you don't want to use our makefile ;) ):

`rivet-buildplugin RivetMYPLUGIN.so MYPLUGIN.cc -lBOOSTFastJets`

If you don't do this you'll still be able to build the analysis,
but Rivet won't be able to access any functions from BOOSTFastJets
you use, and will crash when it attempts to do so. The error message
usually looks something like:

python: symbol lookup error: *bunch of crap*

Physics Motivation
------------------

The big picture aim of this study is to provide an accurate picture of
how different Monte Carlo generators handle the creation of jets, and
in particular the substructure of jets.

Implemented code
------------------------
* Jet Charge
* Jet Dipolarity
* Jet Pull
* N-subjettiness (w/ minimisation algorithm)
* Angular Correlation/Structure Functions
* Grooming (wrapper of FastJet stuff)

Implemented as an example in the SUBSTRUCTURE analysis:

* Average ASF (note this not completely validated, don't rely on it without first making sure you're happy with it!)
