## Makefile to build BOOST2012 MC Rivet analysis

CC=g++
WFLAGS= -Wall -Wextra
CFLAGS=-m64 -pg -I$(INCDIR) -I$(RIVETINCDIR) -O2 $(WFLAGS) -pedantic -ansi

INCDIR=$(PWD)/include
LIBDIR:=$(shell rivet-config --libdir)
PREFIX:=$(shell rivet-config --prefix)
RIVETINCDIR:=$(shell rivet-config --includedir)
LDFLAGS:=$(shell rivet-config --ldflags)

all: rivet-charge rivet-substructure
rivet-charge: libBOOSTFastJets.so
	$(CC) -shared -fPIC $(CFLAGS) -o "RivetMC_GENSTUDY_JETCHARGE.so" MC_GENSTUDY_JETCHARGE.cc -lBOOSTFastJets -L ./ $(LDFLAGS)
rivet-substructure: libBOOSTFastJets.so
	$(CC) -shared -fPIC $(CFLAGS) -o "RivetMC_GENSTUDY_JET_SUBSTRUCTURE.so" MC_GENSTUDY_JET_SUBSTRUCTURE.cc -lBOOSTFastJets -L ./ $(LDFLAGS)
libBOOSTFastJets.so:
	$(CC) -shared -fPIC $(CFLAGS) src/BOOSTFastJets.cxx -o libBOOSTFastJets.so -lfastjet -lfastjettools $(LDFLAGS)
install:
	cp libBOOSTFastJets.so $(LIBDIR)
#	cp RivetMC_GENSTUDY_JETCHARGE.so $(LIBDIR)
#	cp MC_GENSTUDY_JETCHARGE.plot $(PREFIX)/share
#	cp MC_GENSTUDY_JETCHARGE.info $(PREFIX)/share
clean:
	rm -f *.o  *.so
