#
# Makefile
#
#  Makefile for COMPMON
#   mini version of fadc analysis code for Compton Monitorying
#  START UPDATE 1/21/2015 gbf  


#---------------------------------
#  Define location by defining one of the symbols below
LOCATION_JLAB_FRANKLIN = 1
#LOCATION_CMU=1

#
#----------------

# -------------------------------------------------------------------------
#   Sources for all modules
# -------------------------------------------------------------------------

SRC	= mainanal_CompMon.cc fadcTriggered.cc fadcdata.cc vmeauxdata.cc comptonStatus.cc textParams.cc fadcAccums.cc THaEpics.C helicityTracker.cc

OBJ     = $(SRC:.cc=.o) THaCodaFile.o THaCodaData.o

NAME    = compmon

# -------------------------------------------------------------------------
#   These are the standard libraries, include paths and compiler settings
# -------------------------------------------------------------------------
ifdef LOCATION_CMU
BINDIR		= /home/franklin/HallA/newCompton/bin
LIBS		= $(shell $(ROOTSYS)/bin/root-config --libs)
GLIBS		= $(ROOTGLIBS) -L/usr/X11R6/lib -lXpm -lX11
#####EVIO_LIB	= libevio.a
EVIO		= /home/franklin/HallA/newCompton/evio/lib
EVIO_LIB	= $(EVIO)/libevio.so
INCLUDES	= -I /home/franklin/HallA/newCompton/evio/include
endif

ifdef LOCATION_JLAB_FRANKLIN
BINDIR		= /home/compton/franklin/bin
LIBS		= $(shell $(ROOTSYS)/bin/root-config --libs)
GLIBS		= $(ROOTGLIBS) -L/usr/X11R6/lib -lXpm -lX11

#####EVIO_LIB	= libevio.a
EVIO		= /home/compton/evio/Linux-x86_64/lib
EVIO_LIB	= $(EVIO)/libevio.so
endif

ALL_LIBS	= $(EVIO_LIB) $(GLIBS) $(ROOTLIBS)

CC		= g++ -D_FILE_OFFSET_BITS=64


ifdef ONLINE
	ALL_LIBS += $(ONLIBS)
endif

# -------------------------------------------------------------------------
#   Compiler & Linker flags
# -------------------------------------------------------------------------
ROOT_FLAGS      = $(shell $(ROOTSYS)/bin/root-config --cflags)

CXXFLAGS	= -g -Wall $(ROOT_FLAGS) $(LIBS) $(INCLUDES) -DSTANDALONE

# =========================================================================
#   Targets & Rules
# =========================================================================
all: object

object: Makefile $(OBJ)
	$(CC)  $(CXXFLAGS) -o $(BINDIR)/$(NAME) $(OBJ)  $(ALL_LIBS)

clean:
	rm -fv *.o